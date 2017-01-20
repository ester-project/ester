jps=0

; data from the spherical Bessel equation (eq_bessel_full)
; original file name = matrixbandAmsB

nligne,'mat',nl & nl=nl-1
xa=fltarr(nl) & ya=xa
openr,1,'mat'
readf,1,nz3p,ncoef
for k=0,nl-1 do begin
 readf,1,i,j,aa
 xa(k)=i
 ya(k)=j
endfor
close,1

window,0,xsize=800,ysize=800
nx=1 & ny=1
!p.multi=[0,nx,ny]

ang=findgen(16)*(!pi*2/16) & r=0.8
usersym,r*cos(ang),r*sin(ang),/fill

theplot:

if (jps eq 1) then impression,'matrix_pedago',nx,ny,1,/encapsulated

plot,ya,nz3p+1-xa,xrange=[0,nz3p+1],yrange=[0,nz3p+1],xst=5,yst=5,psym=8

; prepare the axis:
tik_val=indgen(nz3p+1)
tik_nam=strarr(nz3p+1)
tik_nam=string(format='(i2)',tik_val)
axis,xaxis=0,xran=[0,nz3p+1],xsty=1,xticks=20,xtickv=tik_val,xtickname=tik_nam
axis,xaxis=1,xran=[0,nz3p+1],xsty=1,xticks=20,xtickv=tik_val,xtickname=tik_nam
axis,yaxis=0,yran=[nz3p+1,0],ysty=1,yticks=20,ytickv=tik_val,ytickname=tik_nam
axis,yaxis=1,yran=[nz3p+1,0],ysty=1,yticks=20,ytickv=tik_val,ytickname=tik_nam
xyouts,12,20,'!94!17Bottom of first domain (bot2)'
xyouts,12,11,'!94!17Top of first domain (top2)'
xyouts,2,11,'!17Top of first domain (top1)!96!17'
xyouts,6,10.2,'(bot1)  !17Bottom of second domain  (bot2)'
xyouts,12,1,'!17Top of second domain (top1)'
xyouts,3,9,'!17Right condition'
xyouts,14,9.2,'!17Left condition'
oplot,[0,40],[10.5,10.5],lin=1
oplot,[10.5,10.5],[0,40],lin=1


if ( jps eq 1 ) then goto,thend

;--------- Postcript ----------------------------------------------
print,'on sort en postscript ? oui=1'
read,oui
if (oui ne 1) then goto,theveryend

jps=1
goto,theplot
;------------------------------------------------------------------

thend:  fin_impression

theveryend: print,'finito !'

end
