!
!     cst_plan95  Elementos finitos elasticos cst = triangulo de
!     ==========  deformacao constante
!
!     Adaptado de:
!     BREBBIA,C.A. & FERRANTE,A.J. (1986) Computational methods for the
!     solution of engineering problems. London: Pentech Press. 370p.
!
!     P¢s-processamento no GiD https://www.gidhome.com/
!
!     Implementado por:
!     Paulo Gustavo Cavalcante Lins <pgcalins@gmail.com>
!
module variaveis
implicit none

  integer, parameter :: ndf=2 !! numero de graus de liberdade por n¢
  integer, parameter :: nne=3 !! numero de n¢s por elemento
  integer, parameter :: ndfel=ndf*nne !! graus de liberdade por elemento

  integer, parameter :: in=15 !! numero do arquivo de entrada
  integer, parameter :: io=16 !! numero do arquivo de saida
  integer, parameter :: i33=33 !! numero do arquivo do GiD
  integer, parameter :: i34=34 !! numero do arquivo do GiD

  integer :: nnode !! numero de n¢s
  real(4), allocatable :: X(:),Y(:) !! coordenadas dos n¢s
  integer :: nelem !! numero de elementos
  integer, allocatable :: kon(:) !! conectividade dos elementos
  integer, allocatable :: imat(:) !! numero do material do elemento
  integer :: nmat !! numero de materiais
  real(4), allocatable :: Ei(:) !! modulo de elasticidade do elemento
  real(4), allocatable :: Poisson(:) !! Poisson do elemento
  real(4), allocatable :: Espessura(:) !! Espessura do elemento
  integer :: nln !! numero de nos carregados
  integer :: nbn !! numero de nos com condicao de contorno
  integer, allocatable :: istatus(:) !! indicador de status
  real(4), allocatable :: Prescrito(:) !! deslocamentos prescritos

  integer :: neq !! numero total de incognitas
  real(4), allocatable :: Carga(:) !! Vetor de cargas nodais
  real(4), allocatable :: Desloc(:) !! Vetor de deslocamentos nodais
  real(4), allocatable :: RigGlobal(:,:) !! Matriz de rigidez global
  real(4) :: RigElem(ndfel,ndfel) !! Matriz de rigidez do elemento

  real(4), allocatable :: Reac(:) !! Tensoes no centroide dos elementos
  real(4), allocatable :: Forc(:) !! Tensoes suavizadas no nos

end module variaveis



module Entrada_Saida

contains

subroutine Abre_Arquivos()
use variaveis

   open(UNIT=in,file='NomeArq.dat',status='old')
   open(UNIT=io,file='NomeArq.out',status='UNKNOWN')
   open(UNIT=i33,file='NomeArq.post.msh',status='UNKNOWN')
   open(UNIT=i34,file='NomeArq.post.res',status='UNKNOWN')

   return
end subroutine Abre_Arquivos



subroutine Entrada_Dados()
use variaveis
implicit none
integer :: i,j,ic(nne),n1,k,k1,L1,L2,n2
real(4) :: w(ndf)

   write(io,'(A,/)')' DADOS DE ENTRADA'

   read(in,*)nnode
   write(io,'(A,i5)')' numeros de nos          :',nnode
   write(io,'(/,A)')' coordenadas nodais'
   write(io,'(7x,A,6x,A,9x,A)')' no ','x','y'
   allocate(X(nnode)); allocate(Y(nnode))
   do j=1,nnode
      read(in,*)i,X(i),Y(i)
      write(io,'(i10,2f10.2)')i,X(i),Y(i)
   enddo

   read(in,*)nelem
   write(io,'(/,A,i5/)')' numero de elementos     :',nelem
   write(io,*)' conectividade dos elementos e propriedades'
   write(io,'(A,3(7x,A),5x,A)')'elemento','no1','no2','no3','material'
   allocate(kon(nne*nelem))
   allocate(imat(nelem))
   do j=1,nelem
      read(in,*)i,ic(1),ic(2),ic(3),imat(i)
      write (io,'(5i10)') i,ic(1),ic(2),ic(3),imat(i)
      n1=nne*(i-1)
      kon(n1+1)=ic(1)
      kon(n1+2)=ic(2)
      kon(n1+3)=ic(3)
   enddo

   read(in,*)nmat
   write(io,'(/,A,i5/)')' numero de materiais     :',nmat
   allocate(Ei(nmat))
   allocate(Poisson(nmat))
   allocate(Espessura(nmat))
   do j=1,nmat
      read(in,*)i,Ei(i),Poisson(i),Espessura(i)
      write(io,'(i10,3f20.5)')i,Ei(i),Poisson(i),Espessura(i)
   enddo

   neq=nnode*ndf  !! calcula o numero total de incognitas
   allocate(Carga(neq))
   allocate(Desloc(neq))
   do i=1,neq
      Carga(i)=0.0 !! Zera vetor de cargas nodais
   enddo

   read(in,*)nln
   write(io,'(/,A,i5/)')' numero de nos carregados:',nln
   write(io,'(A,/7x,A,7x,A,8x,A)')'cargas nodais','no','px','py'
   do i=1,nln
      read(in,*) j,(w(k),k=1,ndf)
      write(io,'(i10,2f10.2)') j,(w(k),k=1,ndf)
      do k=1,ndf
         k1=ndf*(j-1)+k
         Carga(k1)=w(k)
      enddo
   enddo

   read(in,*)nbn
   write(io,'(//,A,i5//)')' numero de nos suportados:',nbn
   write(io,*)' dados das condicoes de contorno'
   write(io,'(23X,A,14X,A)')'status','valores prescritos'
   write(io,'(16X,A)')'(1:prescrito, 0:livre)'
   write(io,'(7x,A,10x,A,9x,A,16x,A,9x,A)')'no','u','v','u','v'
   allocate(istatus((ndf+1)*nbn))
   allocate(Prescrito(neq))
   do i=1,nbn
      read(in,*) j,(ic(k),k=1,ndf),(w(k),k=1,ndf)
      write(io,'(3i10,10x,2f10.4)') j,(ic(k),k=1,ndf),(w(k),k=1,ndf)
      L1=(ndf+1)*(i-1)+1
      L2=ndf*(j-1)
      istatus(L1)=j
      do k=1,ndf
         n1=L1+k
         n2=L2+k
         istatus(n1)=ic(k)
         Prescrito(n2)=w(k)
         !write(io,*)ndf*(nnode-1)+ndf,'L1',L1,'L2',L2,'n1',n1,'n2',n2
      enddo
   enddo

return
end subroutine Entrada_Dados



subroutine Imprime_Resultados()
use variaveis
implicit none
integer :: i,j,k1,k2

      write(io,'(//,A,/)')' RESULTADOS'

      write(io,'(A,/)')' deslocamentos nodais'
      write(io,'(8X,A,12X,A,14x,A)')'no','u','v'
      do i=1,nnode
        k1=ndf*(i-1)+1
        k2=k1+ndf-1
        write(io,'(i10,2E15.5)') i,(Desloc(j),j=k1,k2)
      enddo

      write(io,'(/,A,/)')' tensoes nos elementos'
      write(io,'(2X,A,3(12X,A))')'elemento','s11','s22','s12'
      do i=1,nelem
        k1=3*(i-1)+1
        k2=k1+2
        write(io,'(i10,3E15.5)') i,(forc(j),j=k1,k2)
      enddo

      write(io,'(/,A,/)')' tensoes medias nodais'
      write(io,'(2X,A,3(12X,A))')'node','s11','s22','s12'
      do i=1,nnode
        k1=3*(i-1)+1
        k2=k1+2
        write(io,'(i10,3E15.5)') i,(Reac(j),j=k1,k2)
      enddo

   return
end subroutine Imprime_Resultados



subroutine Resultados_GiD()
use variaveis
implicit none
integer :: i,j,n1,k1,k2

   !open(UNIT=i33,file=arq1,status='UNKNOWN')
   !! Escreve arquivo da malha
   write(i33,*)'MESH "Malha" dimension 2 ElemType Triangle Nnode 3'
   write(i33,*)'Coordinates'
   do i=1,nnode
      write(i33,*)i,'  ',x(i),'  ',y(i)
   enddo
   write(i33,*)'end coordinates'
   write(i33,*)'Elements'
   do i=1,nelem
      n1=nne*(i-1)
      write(i33,'(5(2X,I8))') i,kon(n1+1),kon(n1+2),kon(n1+3),imat(i)
   enddo
   write(i33,*)'end elements'
   !close(UNIT=i33)

   !open(UNIT=i34,file=arq2,status='UNKNOWN')
   !! Escreve arquivo de resultados
   write(i34,*)'GiD Post Results File 1.0'
   write(i34,*)'GaussPoints "Malha_gauss" ElemType Triangle "Malha"'
   write(i34,*)'Number Of Gauss Points: 1'
   write(i34,*)'Nodes not included'
   write(i34,*)'Natural Coordinates: Internal'
   write(i34,*)'End gausspoints'
   write(i34,*)'Result "Desloc" "Load Analysis" ',1,' Vector OnNodes'
   write(i34,*)'ComponentNames "X-Desloc", "Y-Desloc"'
   write(i34,*)'Values'
   do i=1,nnode
      k1=ndf*(i-1)+1
      k2=k1+ndf-1
      write(i34,'(i10,2(2X,E15.5))')i,(Desloc(j),j=k1,k2)
   enddo
   write(i34,*)'End Values'
   write(i34,*)'Result "Gauss Stress" "Load Analysis" ',1,' ',&
              & 'PlainDeformationMatrix OnGaussPoints "Malha_gauss"'
   write(i34,*)'Values'
      do i=1,nelem
        k1=3*(i-1)+1
        k2=k1+2
        write(i34,'(i10,3E15.5)') i,(forc(j),j=k1,k2)
      enddo
   write(i34,*)'End Values'
   write(i34,*)'Result "NodalStress" "Load Analysis" ',1,&
              & ' PlainDeformationMatrix OnNodes "Malha"'
   write(i34,*)'Values'
      do i=1,nnode
        k1=3*(i-1)+1
        k2=k1+2
        write(i34,'(i10,3E15.5)') i,(Reac(j),j=k1,k2)
      enddo
   write(i34,*)'End Values'
   !close(UNIT=i34)

   return
end subroutine Resultados_GiD



subroutine Fecha_Arquivos()
use variaveis

   close(UNIT=in)
   close(UNIT=io)
   close(UNIT=i33)
   close(UNIT=i34)

   return
end subroutine Fecha_Arquivos


end module Entrada_Saida



module Matriz_Elemento

contains

subroutine Monta_Matriz_Elemento(nel)
use variaveis
implicit none
integer :: nel,L,n1,n2,n3
integer :: i,j,IX,IZ,JX,JZ
real(4) :: B(3),D(3),A,ALA,AX,E,G

      L=nne*(nel-1)  !! nel = numero do elemento corrente
      n1=kon(L+1)    !! ponto nodal 1
      n2=kon(L+2)    !! ponto nodal 2
      n3=kon(L+3)    !! ponto nodal 3
      B(1)=Y(N2)-Y(N3)
      B(2)=Y(N3)-Y(N1)
      B(3)=Y(N1)-Y(N2)
      D(1)=X(N3)-X(N2)
      D(2)=X(N1)-X(N3)
      D(3)=X(N2)-X(N1)
      A=(B(1)*D(2)-B(2)*D(1))/2
      do i=1,3
        B(I)=B(I)/(2.*A)
        D(I)=D(I)/(2.*A)
      enddo
      E=Ei(imat(nel))
      G=Poisson(imat(nel))
      ALA=E*Espessura(imat(nel))*A/(1-G**2.)
      AX=(1.-G)/2.
      do i=1,6
        do j=1,6
          RigElem(i,j)=0.0 !! Zera matriz de rigidez do elemento
        enddo
      enddo
      do i=1,3
        IX=(I-1)*2+1
        IZ=(I-1)*2+2
        do J=1,3
          JX=(J-1)*2+1
          JZ=(J-1)*2+2
          RigElem(IX,JX)=(B(I)*B(J)+AX*D(I)*D(J))*ALA
          RigElem(IX,JZ)=(G*B(I)*D(J)+AX*D(I)*B(J))*ALA
          RigElem(IZ,JZ)=(D(I)*D(J)+AX*B(I)*B(J))*ALA
          RigElem(IZ,JX)=(G*D(I)*B(J)+AX*B(I)*D(J))*ALA
        enddo
      enddo

   return
end subroutine Monta_Matriz_Elemento



subroutine Calcula_Incognitas_Secundarias()
use variaveis
implicit none
integer :: i,J,L,nel,n1,n2,n3,K1,K2,K3,L1,L2,L3,L4
real(4) :: B(3),D(3),A,E,ANU,C,ANUP
integer, allocatable :: nodes(:) !! Contador de elementos

      allocate(Reac(3*nnode))
      do i=1,3*nnode
        Reac(i)=0.0
      enddo
      allocate(Forc(3*nelem))
      do i=1,3*nelem
        Forc(i)=0.0
      enddo
      allocate(nodes(nnode))
      do i=1,nnode
        nodes(i)=0
      enddo

      do nel=1,nelem
        L=nne*(nel-1)  !! nel = numero do elemento corrente
        n1=kon(L+1)    !! ponto nodal 1
        n2=kon(L+2)    !! ponto nodal 2
        n3=kon(L+3)    !! ponto nodal 3
        !! Define coordenadas do triangulo
        !! Vetores B e D contem segundo e terceiro elemento da coluna da
        !! matriz C. A = Area do elemento vezes 2
        B(1)=Y(N2)-Y(N3)
        B(2)=Y(N3)-Y(N1)
        B(3)=Y(N1)-Y(N2)
        D(1)=X(N3)-X(N2) !! Primeiro lado do elemento
        D(2)=X(N1)-X(N3) !! Segundo lado do elemento
        D(3)=X(N2)-X(N1) !! Terceiro lado do elemento
        A=(B(1)*D(2)-B(2)*D(1))
        do I=1,3
          B(I)=B(I)/A
          D(I)=D(I)/A
        enddo
        !! Calcula tensoes para o elemento nel
        ANU=Poisson(imat(nel))
        E=Ei(imat(nel))
        ANUP=(1-ANU)/2
        K1=NDF*(N1-1)
        K2=NDF*(N2-1)
        K3=NDF*(N3-1)
        C=E/(1-ANU*ANU)
        L=3*(NEL-1)
        FORC(L+1)=C*(B(1)*Desloc(K1+1)+B(2)*Desloc(K2+1)+B(3)*Desloc(K3+1)+  &
               &  ANU*(D(1)*Desloc(K1+2)+D(2)*Desloc(K2+2)+D(3)*Desloc(K3+2)))
        FORC(L+2)=C*(ANU*(B(1)*Desloc(k1+1)+B(2)*Desloc(K2+1)+B(3)*Desloc(K3+1)) &
              &   +D(1)*Desloc(K1+2)+D(2)*Desloc(K2+2)+D(3)*Desloc(K3+2))
        FORC(L+3)=C*(1-ANU)*(D(1)*Desloc(K1+1)+D(2)*Desloc(K2+1)+D(3)* &
              &  Desloc(K3+1)+B(1)*Desloc(K1+2)+B(2)*Desloc(K2+2)+B(3)*Desloc(K3+2))/2.
        !! Calcula as tensoes nodais medias, no vetor Reac
        K1=3*(N1-1)
        K2=3*(N2-1)
        K3=3*(N3-1)
        do I=1,3
          L1=K1+I
          L2=K2+I
          L3=K3+I
          L4=L+I
          REAC(L1)=REAC(L1)+FORC(L4)
          REAC(L2)=REAC(L2)+FORC(L4)
          REAC(L3)=REAC(L3)+FORC(L4)
        enddo
        !! Vetor nodes contem o numero de elementos concectados com
        !! cada no
        NODES(N1)=NODES(N1)+1
        NODES(N2)=NODES(N2)+1
        NODES(N3)=NODES(N3)+1
      enddo
      do I=1,nnode
        K1=3*I-2
        K2=K1+2
        do J=K1,K2
          REAC(J)=REAC(J)/NODES(I)
        enddo
      enddo

   return
end subroutine Calcula_Incognitas_Secundarias


end module Matriz_Elemento



module Matriz_Global

contains

subroutine Monta_Matriz_Global()
use variaveis
use Matriz_Elemento
implicit none
integer :: i,j,nel

    allocate(RigGlobal(neq,neq))
    do i=1,neq
       do j=1,neq
          RigGlobal(i,j)=0.0  !! Zera matriz global
       enddo
    enddo

    do nel=1,nelem
       call Monta_Matriz_Elemento(nel)
       call Armazena_Elemento_na_Global(nel)
    enddo

   return
end subroutine Monta_Matriz_Global



subroutine Armazena_Elemento_na_Global(nel)
use variaveis
implicit none
integer :: nel
integer :: i,j,i1,j1,i2,j2,n1,n2,k,L,jr,kr,jc,kc

   do i=1,nne
      n1=kon(nne*(nel-1)+i)
      i1=ndf*(i-1)
      j1=ndf*(n1-1)
      do j=1,nne
         n2=kon(nne*(nel-1)+j)
         i2=ndf*(j-1)
         j2=ndf*(n2-1)
         do k=1,ndf
            jr=i1+k
            kr=j1+k
            do L=1,ndf
               jc=i2+L
               kc=j2+L
               RigGlobal(kr,kc)=RigGlobal(kr,kc)+RigElem(jr,jc)
            enddo
         enddo
      enddo
   enddo

   return
end subroutine Armazena_Elemento_na_Global



subroutine Impoe_Condicoes_de_Contorno()
use variaveis
implicit none
integer :: L,no,k1,i,kr,j

   do L=1,nbn
      no=istatus((ndf+1)*(L-1)+1)
      k1=ndf*(no-1)
      do i=1,ndf
         if (istatus((ndf+1)*(L-1)+1+i).eq.1) then
            kr=k1+i
            do j=1,neq
               Carga(j)=Carga(j)-RigGlobal(kr,j)*Prescrito(kr)
               RigGlobal(kr,j)=0.0
               RigGlobal(j,kr)=0.0
            enddo
            RigGlobal(kr,kr)=1.0
            Carga(kr)=Prescrito(kr)
         endif
      enddo
   enddo

   return
end subroutine Impoe_Condicoes_de_Contorno



subroutine Sistema_Linear_Gauss(n,A,B,X,io)
implicit none
integer :: n,io
real(4) :: A(n,n),B(n),X(n)
integer :: i,j,k
real(4) :: c

   do k=1,(n-1)
      c=A(k,k)
      if (abs(c).lt.1.0E-7) then
         write(io,*)'Singularidade na linha ',k
         return
      endif
      do j=1,n
         A(k,j)=A(k,j)/c
      enddo
      B(k)=B(k)/c
      do i=(k+1),n
         c=A(i,k)
         do j=1,n
            A(i,j)=A(i,j)-c*A(k,j)
         enddo
         B(i)=B(i)-c*B(k)
      enddo
   enddo
   if (abs(A(n,n)).lt.1.0E-7) then
      write(io,*)'Singularidade na linha ',k
      return
   endif
   B(n)=B(n)/A(n,n)
   do i=(n-1),1,-1
      do j=(i+1),n
         B(i)=B(i)-A(i,j)*B(j)
      enddo
   enddo
   do i=1,n
      X(i)=B(i)
   enddo

   return
end subroutine Sistema_Linear_Gauss


end module Matriz_Global



program cst_plan95
use variaveis
use Entrada_Saida
use Matriz_Elemento
use Matriz_Global

implicit none

   call Abre_Arquivos()
   call Entrada_Dados()
   call Monta_Matriz_Global()
   call Impoe_Condicoes_de_Contorno()
   call Sistema_Linear_Gauss(neq,RigGlobal,Carga,Desloc,io)
   call Calcula_Incognitas_Secundarias()
   call Imprime_Resultados()
   call Resultados_GiD()
   call Fecha_Arquivos()

end program cst_plan95

