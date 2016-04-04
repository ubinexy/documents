        subroutine legepol_sum(x,n,pol,der,sum)
        implicit real *8 (a-h,o-z)
c
        done=1
        sum=0 
c
        pkm1=1
        pk=x
        sum=sum+pkm1**2 /2
        sum=sum+pk**2 *(1+done/2)
c
        pk=1
        pkp1=x
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200

        sum=0 
c
        pol=1
        der=0
        sum=sum+pol**2 /2
        if(n .eq. 0) return
c
        pol=x
        der=1
        sum=sum+pol**2*(1+done/2)
        return
 1200 continue
c
c       n is greater than 1. conduct recursion
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        sum=sum+pkp1**2*(k+1+done/2)
 2000 continue
c
c        calculate the derivative
c
        pol=pkp1
        der=n*(x*pkp1-pk)/(x**2-1)
        return
        end
