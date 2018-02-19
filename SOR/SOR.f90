subroutine SOR(a, b, x, n, eps, kpp, kp)
    !!! Test Git
    implicit real*8(a-h, o-z)
    dimension a(50, 50), b(50), x(50), xx(50)
    !!! Set optimized value of relaxation factor w
    Pai=3.14
    s=Pai/float(n)
    w=2.0/(1.0+sin(s))
    write(6,5) s,w
5   format(1H, ' s=', 2(F10.5))


    do 10 i=1,n
        xx(i)=0.0
        x(i)=0.0
        !!! Calc. Q
        b(i)=b(i)/a(i,i)
10  continue

    do 11 i=1,n
        t=a(i,i)
        do 12 j=1,n
            !!! Calc. R
            a(i,j)=a(i,j)/t
12      continue
11  continue

    !!!!! Iteration for Solve
    kp=1
25  do 15 i=1,n
        sum=0.0
        do 16 j=1,n
            !!! if(i.eq.j) goto 16
            !!! Calc. 2nd term of Eq.(2-4) or (2-10)
            sum=sum +x(j)*a(i,j)
16      continue
        !!! Calc. 1st term of Eq.(2-4) or (2-10) or (2-26)
        x(i)=b(i)-sum
        x(i)=x(i)*w+xx(i)
15  continue

    !!! Convergence Judgement
    do 21 i=1,n
        if(abs(x(i)-xx(i)).gt.eps) goto 22
21  continue
    goto 35
22  if(kp.ge.kpp) goto 35

    !!! Prepare for next iteration
!!    write(6,6) (x(i),i=1,n)
6   format(1H, 4(5X, F10.5))
    kp=kp+1
    do 23 i=1,n
        xx(i)=x(i)
23  continue
    goto 25

35  return

end subroutine SOR

program Test_SOR
    implicit real*8(a-h,o-z)
    dimension a(50,50),b(50),x(50)
    data (a(1,i),i=1,4) / 7.0, 2.0,-1.0, 1.0/
    data (a(2,i),i=1,4) / 1.0, 5.0, 1.0,-2.0/
    data (a(3,i),i=1,4) / 2.0, 3.0, 8.0, 1.0/
    data (a(4,i),i=1,4) / 2.0,-2.0,-1.0,10.0/
    data (b(i),i=1,4)   /12.0, 6.0,36.0,35.0/
    n=4
    eps=0.000001
    kpp=50

    call SOR(a,b,x,n,eps,kpp,kp)
    write(6,7) kp
7   format(1H, ' P=', I3)
    write(6,8) (x(i),i=1,n)
8   format(1H, 4(5X, F10.5))

    stop
end program Test_SOR

