!-----------------------------!
 subroutine initgraph(frame,t)
!-----------------------------!
 use parameters
 implicit none

 integer :: i1,i2,i3,frame
 real(8) :: t
 character(9) :: fname

 fname='psi000.ps'
 i1=mod(frame,10)
 i2=mod(frame,100)/10
 i3=frame/100
 fname(4:4)=achar(48+i3)
 fname(5:5)=achar(48+i2)
 fname(6:6)=achar(48+i1)
 open(1,file=fname)
 write(1,*)'%!'
 write(1,*)'/m {moveto} def'
 write(1,*)'/l {lineto} def'
 write(1,*)'/TimesRoman findfont 24 scalefont setfont'
 write(1,*)'300 200 translate'
 write(1,*)'-200 0 moveto 400 0 rlineto 0 300 rlineto -400 0 rlineto'
 write(1,*)'closepath stroke'
 write(1,*)'0 0 moveto 0 5 rlineto'
 write(1,*)'-100 0 moveto 0 5 rlineto'
 write(1,*)'100 0 moveto 0 5 rlineto'
 write(1,*)'-200 60 moveto 5 0 rlineto'
 write(1,*)'-200 120 moveto 5 0 rlineto'
 write(1,*)'-200 180 moveto 5 0 rlineto'
 write(1,*)'-200 240 moveto 5 0 rlineto'
 write(1,*)'-50 268 moveto'
 write(1,10)t/fs2au
 write(1,*)'-20 -20 moveto'
 write(1,11) 0.0
 write(1,*)'80 -20 moveto'
 write(1,11) 0.5
 write(1,*)'180 -20 moveto'
 write(1,11) 1.0
 write(1,*)'-120 -20 moveto'
 write(1,11) -0.5
 write(1,*)'-220 -20 moveto'
 write(1,11) -1.0
 write(1,*)'-240 50 moveto'
 write(1,11) 1.0
 write(1,*)'-240 110 moveto'
 write(1,11) 2.0
 write(1,*)'-240 170 moveto'
 write(1,11) 3.0
 write(1,*)'-240 230 moveto'
 write(1,11) 4.0
 write(1,*)'-240 290 moveto'
 write(1,11) 5.0
 write(1,*)'-20 -60 moveto'
 write(1,*) '(x (au)) show'
 write(1,*) 'gsave'
 write(1,*)'-250 130 moveto'
 write(1,*)'((x,t)|) 90 rotate show'
 write(1,*) 'grestore'
 write(1,*) 'gsave'
 write(1,*)'/Symbol findfont 24 scalefont setfont'
 write(1,*)'-250 106 moveto'
 write(1,*)'(|Y) 90 rotate show'
 write(1,*) 'grestore'
 write(1,*) 'gsave'
 write(1,*)'/TimesRoman findfont 16 scalefont setfont'
 write(1,*)'-265 175 moveto'
 write(1,*)'(2) 90 rotate show'
 write(1,*) 'grestore'
 10 format ('(t=',f5.1,' fs) show')
 11 format ('(',f4.1,') show')
 
 end subroutine initgraph 
!------------------------!

!-------------------------------!
 subroutine graphpsi(dx,npoints,psi)
!-------------------------------!
 implicit none

 integer :: i,j,npoints
 real(8) :: x,dx,psi(npoints)
 logical :: firstline

 firstline=.true.
 write(1,*)'gsave 1 0 0 setrgbcolor 2 setlinewidth'
 do i=-npoints/2+1,npoints/2
    if (i>0) then
       j=i
    else     
       j=i+npoints
    endif
    x=dble(i)*dx
    if (abs(x)<=1.d0) then
      if (firstline) then
         write(1,20)200.d0*x,60.d0*psi(j)
         20 format(f7.2,' ',f7.2,' m')
         firstline=.false.
      else
         write(1,10)200.d0*x,60.d0*psi(j)
         10 format(f7.2,' ',f7.2,' l')
      endif
    endif
 end do
 write(1,*)'stroke grestore'
 close(1)

 end subroutine graphpsi
!-----------------------!

!-------------------------------!
 subroutine graphpot(dx,npoints)
!-------------------------------!
 use parameters
 implicit none

 integer :: i,j,npoints
 real(8) :: x,dx,potential,potscale,cutoff
 logical :: firstline

 if (potentialtype=='harmonic') then
    potscale=20.d0
    cutoff=16.d0
 elseif (potentialtype=='doublewell') then
    potscale=40.d0
    cutoff=8.d0
 endif

 firstline=.true.
 write(1,*)'gsave 0 0 1 setrgbcolor 2 setlinewidth'
 do i=-npoints/2+1,npoints/2
    x=dble(i)*dx
    if (potentialtype=='harmonic') then
       potential=0.5d0*mass*angfreq**2*x**2*au2kcalmol
    elseif (potentialtype=='doublewell') then
       potential=barrier*(16.d0*x**4 - 8.d0*x**2 + 1.d0)
    endif
    if (abs(x)<=1.d0.and.potential<cutoff) then
      if (firstline) then
           write(1,20)200.d0*x,potscale*potential
         20 format(f7.2,' ',f7.2,' m')
         firstline=.false.
      else
         write(1,10)200.d0*x,potscale*potential
         10 format(f7.2,' ',f7.2,' l')
      endif
    endif
 end do
 write(1,*)'stroke grestore'

 end subroutine graphpot
!-----------------------!
