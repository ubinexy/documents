        subroutine legewhts(n,ts,whts,ifwhts)
        implicit real *8 (a-h,o-z)
        dimension ts(1),whts(1),ws2(1000),rats(1000)
c
c        this subroutine constructs the nodes and the
c        weights of the n-point gaussian quadrature on 
c        the interval [-1,1]
c
c                input parameters:
c
c  n - the number of nodes in the quadrature
c
c                output parameters:
c
c  ts - the nodes of the n-point gaussian quadrature
c  w - the weights of the n-point gaussian quadrature
c
c       . . . construct the array of initial approximations
c             to the roots of the n-th legendre polynomial
c
        eps=1.0d-14
        ZERO=0
        DONE=1
        pi=datan(done)*4
        h=pi/(2*n) 
        do 1200 i=1,n
        t=(2*i-1)*h
        ts(n-i+1)=dcos(t)
1200  CONTINUE
c
c         use newton to find all roots of the legendre polynomial
c
        ts(n/2+1)=0
        do 2000 i=1,n/2
c
        xk=ts(i)
        ifout=0
        deltold=1
        do 1400 k=1,10
        call legepol_sum(xk,n,pol,der,sum)
        delta=-pol/der
        xk=xk+delta
        if(abs(delta) .lt. eps) ifout=ifout+1
c
        if(ifout .eq. 3) goto 1600
 1400 continue
 1600 continue
        ts(i)=xk
        ts(n-i+1)=-xk
 2000 continue
c     
c        construct the weights via the orthogonality relation
c
        if(ifwhts .eq. 0) return
c
        do 2400 i=1,(n+1)/2
        call legepol_sum(ts(i),n,pol,der,sum)
        whts(i)=1/sum
        whts(n-i+1)=whts(i)
 2400 continue
c
        return
        end
