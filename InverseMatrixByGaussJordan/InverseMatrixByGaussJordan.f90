program InverseMatrixByGaussJordan

    implicit real*8 (a-h,o-z)
    dimension a(10,20)
    data (a(1,i),i=1,3) / 0.0, 2.0, 3.0 /
    data (a(2,i),i=1,3) / 4.0, 0.0, 5.0 /
    data (a(3,i),i=1,3) / 1.0, 2.0, 0.0 /
    n=3
    n1=n+1
    m=n*2
    e=0.000001


    do 10 i=1,n
        do 11 j=n1,m
            a(i,j)=0.0
            if(j.eq.i+n) a(i,j)=1.0
11      continue
10  continue

!    do 999 i=1,n
!        write(6,99) (a(i,j),j=1,n)
!99       format(1H, 5X, 3(3X,E16.8))
!999 continue

! Initialization is Over

    do 15 i=1,n
        ma=i
        fmax=abs(a(i,i))
        if(i.eq.n) goto 50
        i1=i+1

! Pivoting Operation
        do 17 ii=i1,n
            fmax1=abs(a(ii,i))
            if(fmax.ge.fmax1) goto 17
            ma=ii
            fmax=fmax1
17      continue

50      if(fmax.le.e) goto 20
        if(ma.eq.i) goto 51

! Interchange
        do 25 j=i,m
            aa=a(ma,j)
            a(ma,j)=a(i,j)
            a(i,j)=aa
25      continue

! Elementary Operation
51      b=a(i,i)
        do 26 j=i,m
            a(i,j)=a(i,j)/b
26      continue

        do 27 l=1,n
            if(l.eq.i) goto 27
            s=a(l,i)
            do 28 j=1,m
                a(l,j)=a(l,j)-s*a(i,j)
28          continue
27      continue
15  continue

    do 7 i=1,n
        write(6,8) (a(i,j),j=n1,m)
8       format(1H, 5X, 3(3X,E16.8))
7   continue
    goto 35

! Error
20  write(6,9) i,fmax
9   format(1H,' K =',I2,5X,'FMAX=',E14.7)

35  stop
!    end

end program InverseMatrixByGaussJordan
