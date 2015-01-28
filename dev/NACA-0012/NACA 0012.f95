Program NACA12
!Convalidação do Metodo

parameter (n=900,m=600,c=400)
real f(0:8,0:n,0:m),feq(0:8,0:n,0:m),rho(0:n,0:m)
real w(0:8),cx(0:8),cy(0:8),u(0:n,0:m),v(0:n,0:m),image(0:n,0:m),cs, NACA(0:c)
integer i,j,k,x
open(2,file='Naca')

uo=0.1                             !10
obstX = n/3
obstY = m/2
NACA(0)=0
do x=1,c
  NACA(x)=c*(0.6*(0.29690*(x/c)**(1/2)-0.126*(x/c)-0.3516*(x/c)**(2)+0.2843*(x/c)**(3)-0.1015*(x/c)**(4)))
end do
wall=1
do i =0,n
  do j=0,m								!20
    image(i,j)=0.0 
  end do
end do 
do i =obstX,obstX+c                           
  do j =0,m
    if (((j) <= (NACA(i-(n/3))+obstY)).and.((j) >= (obstY-NACA(i-(n/3))))) image(i,j) = wall
  end do
end do                                                                                                                                
rhoo=1.0              
dx=0.25								!30
dy=dx
dt=0.25
cs=((dx)/dt)/1.73205
visco=0.005
Re=uo*c*dx/visco                                        !30
print*,"Re=",Re
omega=1.0/(((3.0*dt*visco)/dx**2)+0.5)
mstep=1500                                       
w(0)=4./9.                                                 
do i=1,4
  w(i)=1./9.
end do
do i=5,8
  w(i)=1./36.
end do                                            !40
cx(0)=0.0
cx(1)=dx/dt
cx(2)=-dx/dt                                             
cx(3)=0.0                                                                 
cx(4)=0.0
cx(5)=dx/dt
cx(6)=-dx/dt
cx(7)=-dx/dt
cx(8)=dx/dt
cy(0)=0.0                                             !50
cy(1)=0.0
cy(2)=0.0
cy(3)=dy/dt                                           
cy(4)=-dy/dt                                                 
cy(5)=dy/dt
cy(6)=-dy/dt
cy(7)=dy/dt
cy(8)=-dy/dt
do j=0,m                                            !60
  do i=0,n
    rho(i,j)=rhoo
    u(i,j)=0.0                                      
    v(i,j)=0.0                                         
  end do
end do
do j=0,m
  do i=0,n
    do k=0,8
      f(k,i,j)=0.0                                !70
    end do
  end do
end do                                               
do j=1,m-1                                             
  u(0,j)=uo
end do


!Loop Principal
do kk=1,mstep                                                 !80
  do i=0,n                     !Processo de Colisão          
    do j=0,m
      t1=u(i,j)*u(i,j)+v(i,j)*v(i,j)                           
      do k=0,8                                                             
        t2=u(i,j)*cx(k)+v(i,j)*cy(k)
        feq(k,i,j)=rho(i,j)*w(k)*(1.0+(t2/cs**(2))+((0.5*t2**2)/cs**4)-(0.5*t1)/cs**2)
        !if(k.eq.0) feq(k,i,j)=w(k)*rho(i,j)
        f(k,i,j)=omega*feq(k,i,j)+(1.0-omega)*f(k,i,j)
      end do
    end do                                                 !90
  end do                                            
  do j=0,m                     !Processo de Transmissão
    do i=1,n                                                       
      f(1,n-i+1,j)=f(1,n-i,j)        !Left to Right             
      f(2,i-1,j)=f(2,i,j)            !Right to Left
    end do                          
  end do
  do i=0,n
    do j=1,m
      f(3,i,m-j+1)=f(3,i,m-j)        !Bottom to Top             100
      f(4,i,j-1)=f(4,i,j)            !Top to Bottom         
    end do
  end do                                                          
  do j=1,m                                                     
    do i=1,n
      f(5,n-i+1,m-j+1)=f(5,n-i,m-j)   !Diagonal ++     
      f(6,i-1,j-1)=f(6,i,j)           !Diagonal --
    end do
  end do
  do j=1,m                                                     !110
    do i=1,n                                               
      f(7,i-1,m-j+1)=f(7,i,m-j)       !Diagonal -+
      f(8,n-i+1,j-1)=f(8,n-i,j)       !Diagonal +-                   
    end do                                                     
  end do
  do  j=0,m                 !Condições de Contorno
    f(1,n,j)=2.*f(1,n-1,j)-f(1,n-2,j)        !Metodo Bounce Back na Camada Leste     
    f(5,n,j)=2.*f(5,n-1,j)-f(5,n-2,j)
    f(8,n,j)=2.*f(8,n-1,j)-f(8,n-2,j)                                                                                      
  end do                                                                    !120
  do i=0,n
    f(3,i,m)=2.*f(3,i,m-1)-f(3,i,m-2)        !Metodo Bounce Back na Camada Norte
    f(5,i,m)=2.*f(5,i,m-1)-f(5,i,m-2)
    f(7,i,m)=2.*f(7,i,m-1)-f(7,i,m-2)                                     
    f(4,i,0)=2.*f(4,i,1)-f(4,i,2)        !Metodo Bounce Back na Camada Sul      
    f(6,i,0)=2.*f(6,i,1)-f(6,i,2)
    f(8,i,0)=2.*f(8,i,1)-f(8,i,2)
  end do
  do j=1,m-1                     !Metodo Bounce Back na Camada Oeste
    rhon=(f(0,0,j)+f(3,0,j)+f(4,0,j)+2.*(f(2,0,j)+f(7,0,j)+f(6,0,j)))/(1.-uo)          !130
    f(1,0,j)=f(2,0,j)+(2.*rhon*uo/3)
    f(5,0,j)=f(6,0,j)+(rhon*uo/6.0)                                   
    f(8,0,j)=f(7,0,j)+(rhon*uo/6.0)                                                
  end do                                                                   
  do i =0,n
    do j =0,m
      if (image(i,j) == wall) then
        if ((image(i,j+1) == wall).and.(image(i,j-1) == wall).and.(image(i-1,j) == wall).and.(image(i+1,j) == 0.0)) then
          f(1,i,j)=f(2,i,j)
          f(5,i,j)=f(6,i,j)       !Direita pra esquerda                                  140
          f(8,i,j)=f(7,i,j)
        end if
        if ((image(i,j+1) == wall).and.(image(i,j-1) == wall).and.(image(i+1,j) == wall).and.(image(i-1,j) == 0.0)) then
          f(2,i,j)=f(1,i,j)
          f(6,i,j)=f(5,i,j)        !Esquerda pra direita
          f(7,i,j)=f(8,i,j)
        end if
        if ((image(i,j+1) == 0.0).and.(image(i,j-1) == wall).and.(image(i-1,j) == wall).and.(image(i+1,j) == wall)) then
          f(3,i,j)=f(4,i,j)
          f(5,i,j)=f(6,i,j)                 !Cima pra baixo                                150
          f(7,i,j)=f(8,i,j)
        end if
        if ((image(i,j+1) == wall).and.(image(i,j-1) == 0.0).and.(image(i-1,j) == wall).and.(image(i+1,j) == wall)) then
          f(4,i,j)=f(3,i,j)
          f(6,i,j)=f(5,i,j)                 !Baixo pra cima
          f(8,i,j)=f(7,i,j)
        end if
        if ((image(i,j+1) == 0.0).and.(image(i,j-1) == wall).and.(image(i-1,j) == wall).and.(image(i+1,j) == 0.0)) then
          f(3,i,j)=f(4,i,j)
          f(5,i,j)=f(6,i,j)                 !Quina direita pra esquerda e cima pra baixo-|              160
          f(7,i,j)=f(8,i,j)
          f(1,i,j)=f(2,i,j)
        end if
        if ((image(i+1,j+1) == 0.0).and.(image(i-1,j-1) == wall).and.(image(i-1,j+1) == wall).and.(image(i+1,j-1) == wall))then
          f(3,i,j)=f(4,i,j)
          f(5,i,j)=f(6,i,j)                 !Quina Inversa direita pra esquerda e cima pra baixo              
          f(7,i,j)=f(8,i,j)
          f(1,i,j)=f(2,i,j)
        end if
        if ((image(i,j+1) == wall).and.(image(i,j-1) == 0.0).and.(image(i-1,j) == wall).and.(image(i+1,j) == 0.0)) then !170
          f(4,i,j)=f(3,i,j)
          f(6,i,j)=f(5,i,j)                 !Quina direita pra esquerda e baixo pra cima_|
          f(8,i,j)=f(7,i,j)
          f(1,i,j)=f(2,i,j)
        end if
        if ((image(i+1,j+1) == wall).and.(image(i-1,j-1) == wall).and.(image(i-1,j+1) == wall).and.(image(i+1,j-1) == 0.0)) then
          f(4,i,j)=f(3,i,j)
          f(6,i,j)=f(5,i,j)                 !Quina Inversa direita pra esquerda e baixo pra cima
          f(8,i,j)=f(7,i,j)
          f(1,i,j)=f(2,i,j)                                                                 !180
        end if
        if ((image(i,j+1) == 0.0).and.(image(i,j-1) == wall).and.(image(i-1,j) == 0.0).and.(image(i+1,j) == wall)) then   
          f(2,i,j)=f(1,i,j)
          f(6,i,j)=f(5,i,j)                 !Quina esquerda pra direita e cima pra baixo |-
          f(7,i,j)=f(8,i,j)
          f(3,i,j)=f(4,i,j)
        end if
        if ((image(i+1,j+1) == wall).and.(image(i-1,j-1) == wall).and.(image(i-1,j+1) == 0.0).and.(image(i+1,j-1) == wall)) then
          f(2,i,j)=f(1,i,j)
          f(6,i,j)=f(5,i,j)                 !Quina Inversa esquerda pra direita e cima pra baixo              190
          f(7,i,j)=f(8,i,j)
          f(3,i,j)=f(4,i,j)
        end if
        if ((image(i,j+1) == wall).and.(image(i,j-1) == 0.0).and.(image(i-1,j) == 0.0).and.(image(i+1,j) == wall)) then
          f(2,i,j)=f(1,i,j)
          f(6,i,j)=f(5,i,j)                 !Quina esquerda pra direita e baixo pra cima|_
          f(7,i,j)=f(8,i,j)
          f(4,i,j)=f(3,i,j)                                                           
        end if
        if ((image(i+1,j+1) == wall).and.(image(i-1,j-1) == 0.0).and.(image(i-1,j+1) == wall).and.(image(i+1,j-1) == wall)) then !200
          f(2,i,j)=f(1,i,j)
          f(6,i,j)=f(5,i,j)                 !Quina esquerda pra direita e baixo pra cima
          f(7,i,j)=f(8,i,j)
          f(4,i,j)=f(3,i,j)                                                           
        end if
      end if                                                           
    end do
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
    do j=0,m
      usum=0.0
      vsum=0.0
      do k=0,8                                                         !160
        usum=usum+f(k,i,j)*cx(k)
        vsum=vsum+f(k,i,j)*cy(k)                                            
      end do
      u(i,j)=usum/rho(i,j)
      v(i,j)=vsum/rho(i,j)
    end do
  end do
  do j=1,m-1
    v(n,j)=0.0
  end do                                                                !170
  do i=0,n
    v(i,m)=0.0
    v(i,0)=0.0
  end do
  do i =0,n
    do j =0,m
      if (image(i,j) == wall) then
      u(i,j)=0.0
      v(i,j)=0.0  
      end if                                                           
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
end program NACA12