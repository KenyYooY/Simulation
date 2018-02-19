program GaussJordan

!   ### Solution of Simultaneous Equation by Gauss-Jordan Method

    implicit real*8 (a-h,o-z)
    dimension a(10,11),x(10)

    data (a(1,i),i=1,5) / 0.0 , 1.0 , 3.0 ,-5.0 ,-9.0 /
    data (a(2,i),i=1,5) / 3.0 ,-2.0 , 1.0 , 1.0 , 6.0 /
    data (a(3,i),i=1,5) / 8.0 ,-5.0 ,-7.0 , 2.0 ,-15.0 /
    data (a(4,i),i=1,5) / 1.0 , 1.0 ,-2.0 ,-3.0 ,-15.0 /

    n=4
    n1=n+1
    e=0.000001

    do 15 i=1,n
        ma=i
        fmax=abs(a(i,i))
        ! print *, fmax
        if(i.eq.n) goto 50
        i1=i+1

        do 17 ii=i1,n
            fmax1=abs(a(ii,i))
            if(fmax.ge.fmax1) goto 17
            ma=ii
            fmax=fmax1
17      continue

50      if(fmax.le.e) goto 20
        if(ma.eq.i) goto 51

!   ### Interchange
        do 25 j=i,n1
            aa=a(ma,j)
            a(ma,j)=a(i,j)
            a(i,j)=aa
25      continue

51      b=a(i,i)

        do 26 j=1,n1
            a(i,j)=a(i,j)/b
26      continue

        do 27 l=1,n
            if(l.eq.i) goto 27
            s=a(l,i)
            do 28 j=1,n1
                a(l,j)=a(l,j)-s*a(i,j)
28          continue
27      continue

15  continue

    do 80 i=1,n
        x(i)=a(i,n1)
80  continue

    write(6,100) (i,x(i),i=1,n)
100 format(//1H ,4(2X,'X',i1,'=',E14.7)//)
    goto 35

20  write(6,200) i,fmax
200 format(1H ,'K=',I2,5X,'FMAX=',F14.7)

35  stop

!    end

end program GaussJordan

