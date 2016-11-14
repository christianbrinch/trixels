      SUBROUTINE fourier_o(nc,nt,nd,X,Y,Z,TR)
c From Lee and Mittra, 1983
c THIS SUBROUTINE WORKS! DONT FIX ANYTHING. CB, 15/3/11
      INTEGER nc,nt,nd
      INTEGER TR(3*nt)
      REAL*4 X(nc), Y(nc), Z(nc)
      COMPLEX S(nd,nd),iu
      PARAMETER (iu=(0,1))
      PARAMETER (pi=3.141592653589793238)
      INTEGER i,j,t,n,m
      DOUBLE PRECISION g(4,3),u,v,p(3),a(3,3),b(3,3),mean,vector(3),w(3)

c Initialize array with zeros
      DO i=1,nd
        DO j=1,nd
          S(i,j)=(0.,0.)
        END DO
      END DO

c loop over the trixels
      DO t=1,nt
        mean=0.

c put triangle vertices in g and calculate mean intensity
        DO n=1,3
          g(n,1)=X(TR(n+3*(t-1))+1)
          g(n,2)=Y(TR(n+3*(t-1))+1)
          g(n,3)=Z(TR(n+3*(t-1))+1)
        END DO
        g(4,1:3)=g(1,1:3)
        mean=(g(1,3)+g(2,3)+g(3,3))/3.

c calculate the slope p of triangle edges
        DO n=1,3
          p(n)=(g(n+1,2)-g(n,2))/(g(n+1,1)-g(n,1))
          vector=g(n+1,1:3)-g(n,1:3)
          a(n,1:3)=vector/DSQRT(DOT_PRODUCT(vector,vector))
          b(n,1:3)=(/ -a(n,2), a(n,1), 0d0 /)
        END DO

c loop over fourier components (u,v)
        DO i=1,nd
          DO j=1,nd
            u=(i-1-nd/2.)*2.*pi*1./nd
            v=(j-1-nd/2.)*2.*pi*1./nd
            w=(/u,v,0d0/)

c loop over triangle vertices
            DO n=1,3
              m=n-1
              IF(m.lt.1) THEN
                m=3
              END IF
c calculate Fourier integral
c              S(i,j)=S(i,j)+
c     &          mean*CDEXP(iu*(u*g(n,1)+v*g(n,2)))*
c     &          (p(m)-p(n))/((u+p(m)*v)*(u+p(n)*v))

              S(i,j)=S(i,j)+
     &          mean*CDEXP(iu*DOT_PRODUCT(w,g(n,1:3)))*
     &          DOT_PRODUCT(b(n,1:3),a(m,1:3))/(DOT_PRODUCT(w,a(n,1:3))
     &          *DOT_PRODUCT(w,a(m,1:3)))

            END DO
          END DO
        END DO
      END DO

c output Fourier components
      DO i=1,nd
        DO j=1,nd
          write(*,*) S(i,j)
        END DO
      END DO

      RETURN
      END






      SUBROUTINE fourier_housh(nc,nt,nd,X,Y,Z,TR)
c From Housmand, Chew, and Lee, 1991
      INTEGER nc,nt,nd
      INTEGER TR(3*nt)
      REAL*4 X(nc), Y(nc), Z(nc)
      COMPLEX S(nd,nd),Vis(nd,nd),f1,f2,gn(3),iu
      PARAMETER (iu=(0,1))
      PARAMETER (pi=3.141592653589793238)
      INTEGER i,j,t,n,m,o
      DOUBLE PRECISION g(4,3),u,v,p(3),a(3,3),b(3,3),mean
      DOUBLE PRECISION vector(3),w(3),h(3)


c Initialize array with zeros
      DO i=1,nd
        DO j=1,nd
          S(i,j)=(0.,0.)
          Vis(i,j)=(0.,0.)
        END DO
      END DO

c loop over the trixels
      DO t=1,nt
        mean=0.

c put triangle vertices in g and calculate mean intensity
        DO n=1,3
          g(n,1)=X(TR(n+3*(t-1))+1)
          g(n,2)=Y(TR(n+3*(t-1))+1)
          g(n,3)=Z(TR(n+3*(t-1))+1)
        END DO
        g(4,1:3)=g(1,1:3)
        mean=(g(1,3)+g(2,3)+g(3,3))/3.

        DO n=1,3
          m=n+1
          IF(m.gt.3) THEN
            m=1
          END IF
          vector=g(m,1:3)-g(n,1:3)
          a(n,1:3)=vector/DSQRT(DOT_PRODUCT(vector,vector))
          b(n,1:3)=(/ -a(n,2), a(n,1), 0d0 /)
        END DO

        DO n=1,3
          m=n+1
          IF(m.gt.3) THEN
            m=1
          END IF
          o=n-1
          IF(o.lt.1) THEN
            o=3
          END IF
          h(n)=DSQRT((g(o,1)-g(n,1))**2+(g(o,2)-g(n,2))**2)
     &         *SIN(ACOS(DOT_PRODUCT(-a(m,1:3),a(o,1:3))))
        END DO


c loop over fourier components (u,v)
        DO i=1,nd
          DO j=1,nd
            u=(i-1-nd/2.)*2.*pi*1./nd
            v=(j-1-nd/2.)*2.*pi*1./nd
            w=(/u,v,0d0/)

c calculate Fourier integral of a flat triangle
            DO n=1,3
              m=n-1
              IF(m.lt.1) THEN
                m=3
              END IF
              S(i,j)=S(i,j)+
     &          mean*CDEXP(iu*DOT_PRODUCT(w,g(n,1:3)))*
     &          DOT_PRODUCT(b(n,1:3),a(m,1:3))/(DOT_PRODUCT(w,a(n,1:3))
     &          *DOT_PRODUCT(w,a(m,1:3)))
            END DO ! end loop over n

            DO n=1,3
              m=n-1
              IF(m.lt.1) THEN
                m=3
              END IF
              o=n+1
              IF(o.gt.3) THEN
                o=1
              END IF

              f1=DOT_PRODUCT(a(o,1:3),b(m,1:3))*
     &           DOT_PRODUCT(a(o,1:3),a(m,1:3))/
     &           (DOT_PRODUCT(w,a(m,1:3)))*
     &           ((1+iu*DOT_PRODUCT(w,(g(n,1:3)-g(m,1:3))))*
     &           CDEXP(-iu*DOT_PRODUCT(w,g(n,1:3)))-
     &           CDEXP(-iu*DOT_PRODUCT(w,g(m,1:3))))

              f2=DOT_PRODUCT(a(o,1:3),b(n,1:3))*
     &           DOT_PRODUCT(a(o,1:3),a(n,1:3))/
     &           (DOT_PRODUCT(w,a(n,1:3)))*
     &           ((1+iu*DOT_PRODUCT(w,(g(n,1:3)-g(o,1:3))))*
     &           CDEXP(-iu*DOT_PRODUCT(w,g(n,1:3)))-
     &           CDEXP(-iu*DOT_PRODUCT(w,g(o,1:3))))

              IF(ABS(DOT_PRODUCT(w,b(o,1:3))).gt.1e-4) THEN
                gn(n)= (-1)/(iu*h(n)*DOT_PRODUCT(w,b(o,1:3)))*
     &                 (S(i,j)+f1+f2)
              END IF
            END DO ! end loop over n

            DO n=1,3
c              Vis(i,j)=Vis(i,j)+g(n,3)*gn(n)
c     &                *CDEXP(iu*DOT_PRODUCT(w,g(n,1:3)))
              Vis(i,j)=Vis(i,j)+g(n,3)*gn(n)
     &                *CDEXP(iu*DOT_PRODUCT(w,g(n,1:3)))
           END DO ! end loop ever n
          END DO ! end loop over j
        END DO ! end loop over i
      END DO ! end loop over t



c output Fourier components
      DO i=1,nd
        DO j=1,nd
          write(*,*) Vis(i,j)
        END DO
      END DO

      RETURN
      END






      SUBROUTINE fourier(nc,nt,nd,VX,VY,VZ,TR)
c From McInturff and Simon, 1991
      INTEGER nc,nt,nd,os
      INTEGER TR(3*nt)
      REAL*4 VX(nc), VY(nc), VZ(nc)
      COMPLEX Vis(nd,nd),iu, T(3,3)
      PARAMETER (iu=(0,1))
      PARAMETER (pi=3.141592653589793238)
      INTEGER i,j,n,tc
      DOUBLE PRECISION g(4,3),u,v,a(3,3),area,w(3),z(3),rnc(3,3)
      DOUBLE PRECISION dp0,dp1,dp2,dp3,dp4,dp5,dp6
      PARAMETER (z=(/0,0,1/))


! Initialize array with zeros
      DO i=1,nd
        DO j=1,nd
          Vis(i,j)=(0.,0.)
        END DO
      END DO

! loop over the trixels
      DO tc=1,nt

! put triangle vertices in g and calculate mean intensity
        DO n=1,3
          g(n,1)=VX(TR(n+3*(tc-1))+1)
          g(n,2)=VY(TR(n+3*(tc-1))+1)
          g(n,3)=VZ(TR(n+3*(tc-1))+1)
        END DO
        g(4,1:3)=g(1,1:3)

        DO n=1,3
          rnc(n,1:3)=(g(n,1:3)+g(n+1,1:3))/2.
          a(n,1:3)=(/ g(n+1,1)-g(n,1), g(n+1,2)-g(n,2), 0.d0 /)
        END DO

        area= abs((g(1,1)-g(2,1))*(g(1,2)-g(3,2))
     &           -(g(1,2)-g(2,2))*(g(1,1)-g(3,1)))


! loop over fourier components (u,v)
        DO i=1,nd
          DO j=1,nd
            u=(i-0.5-nd/2.)*2.*pi*1./nd
            v=(j-0.5-nd/2.)*2.*pi*1./nd


            IF(abs(u).lt.1d-10.and.abs(v).lt.1d-10) THEN
              u = 1d-10
              v = 1d-10
            ENDIF
            w=(/u,v,0d0/)


            DO n=1,3
              T(n,1:3)=(/ (0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0)/)
              DO m=1,3
                dp0=w(1)*rnc(m,1) + w(2)*rnc(m,2) + w(3)*rnc(m,3)
                dp1=w(1)*w(1)+w(2)*w(2)+w(3)*w(3)
                dp2=z(1)*(-w(2)*a(m,3)) + z(2)*w(1)*a(m,3) +
     &           z(3)*a(m,1)*w(2)-w(1)*a(m,2)
                dp3=w(1)*a(m,1) + w(2)*a(m,2) + w(3)*a(m,3)
                dp4=z(1)*(-w(2)*a(m,3)) + z(2)*w(1)*a(m,3) + z(3) +
     &           a(m,1)*w(2)-w(1)*a(m,2)
                dp5=w(1)*a(m,1) + w(2)*a(m,2) + w(3)*a(m,3)

              T(n,1:3)=T(n,1:3)+1./dp1*CDEXP(iu*dp0)*(((/-a(m,2),a(m,1),
     &         0.d0/)+(iu*rnc(m,1:3)-iu*g(n+1,1:3)-2.*w/dp1)*dp2)*
     &         BesJ0(dp3/2)-a(m,1:3)*dp4/2*BesJ1(dp5/2))


c              T(n,1:3)=T(n,1:3)+1./DOT_PRODUCT(w,w)*
c     &                   CDEXP(iu*DOT_PRODUCT(w,rnc(m,1:3))) *
c     &                   (((/-a(m,2), a(m,1), 0.d0/)+(iu*rnc(m,1:3)-
c     &                   iu*g(n+1,1:3)-2.*w/DOT_PRODUCT(w,w)) *
c     &                   DOT_PRODUCT(z,(/-w(2)*a(m,3),w(1)*a(m,3),
c     &                   a(m,1)*w(2)-w(1)*a(m,2)/)))
c     &                   *BesJ0(DOT_PRODUCT(w,a(m,1:3))/2.)
c     &                   -a(m,1:3)*DOT_PRODUCT(z,(/-w(2)*a(m,3),w(1)*
c     &                   a(m,3),a(m,1)*w(2)-w(1)*a(m,2)/) )/2.
c     &                   *BesJ1(DOT_PRODUCT(w,a(m,1:3))/2.))
              END DO ! end loop over m
            END DO ! end loop over n

c          dp1=-a(2,2)*T(2,1) + a(2,1)*T(2,2)
c          dp2=-a(3,2)*T(3,1) + a(3,1)*T(3,2)
c          dp3=-a(1,2)*T(1,1) + a(1,1)*T(1,2)
c
c          VIS(i,j)=VIS(i,j)+(g(1,3)*dp1+g(2,3)*dp2+g(3,3)*dp3
c     &                      )/(2.*area) * ((3.14159/2.)*nd)**2

          VIS(i,j) = VIS(i,j)+
     &             (
     & g(1,3)*DOT_PRODUCT((/-a(2,2), a(2,1),0.d0/)/(2.*area),T(2,1:3)) +
     & g(2,3)*DOT_PRODUCT((/-a(3,2), a(3,1),0.d0/)/(2.*area),T(3,1:3)) +
     & g(3,3)*DOT_PRODUCT((/-a(1,2), a(1,1),0.d0/)/(2.*area),T(1,1:3))
     &             ) !* ((3.14159/2.)*nd)**2


          END DO
        END DO
      END DO



      OPEN (7, FILE = 'temp_in', ACCESS = 'APPEND',STATUS = 'NEW')
      OPEN (8, FILE = 'temp_co', ACCESS = 'APPEND',STATUS = 'NEW')
! output Fourier components
      DO i=1,nd
        DO j=1,nd
           u=(i-0.5-nd/2.)*2.*3.14159*1./nd
           v=(j-0.5-nd/2.)*2.*3.14159*1./nd
           write(8,*) u,v
           write(7,*) REAL(Vis(i,j)), AIMAG(Vis(i,j))
        END DO
      END DO

      CLOSE(7)
      CLOSE(8)
      RETURN
      END
