      subroutine ctang(ndp,xd,yd,nout,nt,ipt,ipl,iwl,iwp,wk)
      implicit double precision (a-h,o-z)
c it carries out triangulation by dividing the x-y plane into a
c number of triangles according to the given data points in the plane.
c at the end of the operations,the indices of the vertexes of the
c triangles are listed counter-clockwise.
c subroutine ctang calls the maxmn function.
c
c
c the input parameters are
c     ndp = number of data points,
c     xd  = array of dimension ndp containing the x coordinates of the
c           data points,
c     yd  = array of dimension ndp containing the y coordinates of the
c           data points,
c     nout= logical unit number for the standard output unit of the
c           system.
c
c the output parameters are
c     nt  = number of triangles,
c     ipt = integer array of dimension 6*ndp-15, where the indices of
c           the vertexes of the (it)th triangle are to be stored as the
c           (3*it-2)nd, (3*it-1)st and (3*it)th elements, it=1,2,..,nt.
c
c the other parameters are
c     ipl = integer array of dimension 6*ndp used internally as work
c           area,
c     iwl = integer array of dimension 18*ndp used internally as work
c           area,
c     iwp = integer array of dimension ndp used internally as work
c           area,
c     wk  = array of dimension ndp used internally as work area.
c
c declaration statements.
      dimension xd(*),yd(*),ipt(*),ipl(*),iwl(*),iwp(*),wk(*)
      dimension itf(2)
      data  ratio/1.0d-6/, nrep/100/
c statement functions.
      dsqf(u1,v1,u2,v2)=(u2-u1)**2+(v2-v1)**2
      side(u1,v1,u2,v2,u3,v3)=(v3-v1)*(u2-u1)-(u3-u1)*(v2-v1)
c preliminary processing.
      ndp0=ndp
      ndpm1=ndp0-1
      if(ndp0.lt.4)       go to 90
c determines the closest pair of data points and their midpoints.
      dsqmn=dsqf(xd(1),yd(1),xd(2),yd(2))
      ipmn1=1
      ipmn2=2
      do 22  ip1=1,ndpm1
        x1=xd(ip1)
        y1=yd(ip1)
        ip1p1=ip1+1
        do 21  ip2=ip1p1,ndp0
          dsqi=dsqf(x1,y1,xd(ip2),yd(ip2))
c         if(dsqi.eq.0.0d0)go to 91
          if(abs(dsqi-0.0d0).le.1.d-10)go to 91
          if(dsqi.ge.dsqmn)go to 21
          dsqmn=dsqi
          ipmn1=ip1
          ipmn2=ip2
   21   continue
   22 continue
      dsq12=dsqmn
      xdmp=(xd(ipmn1)+xd(ipmn2))/2.0d0
      ydmp=(yd(ipmn1)+yd(ipmn2))/2.0d0
c sorts the other (ndp-2) data points in ascending order of distance
c from the midpoints and stores the stored data points numbers in the
c iwp array.
      jp1=2
      do 31  ip1=1,ndp0
        if(ip1.eq.ipmn1.or.ip1.eq.ipmn2)      go to 31
        jp1=jp1+1
        iwp(jp1)=ip1
        wk(jp1)=dsqf(xdmp,ydmp,xd(ip1),yd(ip1))
   31 continue
      do 33  jp1=3,ndpm1
        dsqmn=wk(jp1)
        jpmn=jp1
        do 32  jp2=jp1,ndp0
          if(wk(jp2).ge.dsqmn)      go to 32
          dsqmn=wk(jp2)
          jpmn=jp2
   32   continue
        its=iwp(jp1)
        iwp(jp1)=iwp(jpmn)
        iwp(jpmn)=its
        wk(jpmn)=wk(jp1)
   33 continue
c if necessary, modifies the ordering so that the
c first three data points are not collinear.
      ar=dsq12*ratio
      x1=xd(ipmn1)
      y1=yd(ipmn1)
      dx21=xd(ipmn2)-x1
      dy21=yd(ipmn2)-y1
      do 36  jp=3,ndp0
        ip=iwp(jp)
        if(dabs((yd(ip)-y1)*dx21-(xd(ip)-x1)*dy21).gt.ar)
     1               go to 37
   36 continue
      go to 92
   37 if(jp.eq.3)    go to 40
      jpmx=jp
      jp=jpmx+1
      do 38  jpc=4,jpmx
        jp=jp-1
        iwp(jp)=iwp(jp-1)
   38 continue
      iwp(3)=ip
c forms the first triangle, stores point numbers of the vertexes of
c the triangles in the ipt array, and stores point numbers of the
c border line segments and the triangle number in the ipl array.
   40 ip1=ipmn1
      ip2=ipmn2
      ip3=iwp(3)
      if(side(xd(ip1),yd(ip1),xd(ip2),yd(ip2),xd(ip3),yd(ip3))
     1     .ge.0.0d0)       go to 41
      ip1=ipmn2
      ip2=ipmn1
   41 nt0=1
      ntt3=3
      ipt(1)=ip1
      ipt(2)=ip2
      ipt(3)=ip3
      nl0=3
      nlt3=9
      ipl(1)=ip1
      ipl(2)=ip2
      ipl(3)=1
      ipl(4)=ip2
      ipl(5)=ip3
      ipl(6)=1
      ipl(7)=ip3
      ipl(8)=ip1
      ipl(9)=1
c adds the remaining (ndp-3) data points, one by one.
      do 79  jp1=4,ndp0
        ip1=iwp(jp1)
        x1=xd(ip1)
        y1=yd(ip1)
c determines the visible line segments.
        ip2=ipl(1)
        jpmn=1
        dxmn=xd(ip2)-x1
        dymn=yd(ip2)-y1
        dsqmn=dxmn**2+dymn**2
        armn=dsqmn*ratio
        jpmx=1
        dxmx=dxmn
        dymx=dymn
        dsqmx=dsqmn
        armx=armn
        do 52  jp2=2,nl0
          ip2=ipl(3*jp2-2)
          dx=xd(ip2)-x1
          dy=yd(ip2)-y1
          ar=dy*dxmn-dx*dymn
          if(ar.gt.armn)       go to 51
          dsqi=dx**2+dy**2
          if(ar.ge.(-armn).and.dsqi.ge.dsqmn)      go to 51
          jpmn=jp2
          dxmn=dx
          dymn=dy
          dsqmn=dsqi
          armn=dsqmn*ratio
   51     ar=dy*dxmx-dx*dymx
          if(ar.lt.(-armx))    go to 52
          dsqi=dx**2+dy**2
          if(ar.le.armx.and.dsqi.ge.dsqmx)    go to 52
          jpmx=jp2
          dxmx=dx
          dymx=dy
          dsqmx=dsqi
          armx=dsqmx*ratio
   52   continue
        if(jpmx.lt.jpmn)  jpmx=jpmx+nl0
        nsh=jpmn-1
        if(nsh.le.0)      go to 60
c shifts (rotates) the ipl array so that the invisible border line
c segments are contained in the first part of the ipl array.
        nsht3=nsh*3
        do 53  jp2t3=3,nsht3,3
          jp3t3=jp2t3+nlt3
          ipl(jp3t3-2)=ipl(jp2t3-2)
          ipl(jp3t3-1)=ipl(jp2t3-1)
          ipl(jp3t3)  =ipl(jp2t3)
   53   continue
        do 54  jp2t3=3,nlt3,3
          jp3t3=jp2t3+nsht3
          ipl(jp2t3-2)=ipl(jp3t3-2)
          ipl(jp2t3-1)=ipl(jp3t3-1)
          ipl(jp2t3)  =ipl(jp3t3)
   54   continue
        jpmx=jpmx-nsh
c adds triangles to the ipt array, updates border line segments in
c the ipl array, and sets flags for the border line segments to be
c reexamined in the iwl array.
   60   jwl=0
        do 64  jp2=jpmx,nl0
          jp2t3=jp2*3
          ipl1=ipl(jp2t3-2)
          ipl2=ipl(jp2t3-1)
          it  =ipl(jp2t3)
c adds a triangle to the ipt array.
          nt0=nt0+1
          ntt3=ntt3+3
          ipt(ntt3-2)=ipl2
          ipt(ntt3-1)=ipl1
          ipt(ntt3)  =ip1
c updates the border line segments in the ipl array.
          if(jp2.ne.jpmx)      go to 61
          ipl(jp2t3-1)=ip1
          ipl(jp2t3)  =nt0
   61     if(jp2.ne.nl0)       go to 62
          nln=jpmx+1
          nlnt3=nln*3
          ipl(nlnt3-2)=ip1
          ipl(nlnt3-1)=ipl(1)
          ipl(nlnt3)  =nt0
c determines the vertex that does not lie on the border line segments.
   62     itt3=it*3
          ipti=ipt(itt3-2)
          if(ipti.ne.ipl1.and.ipti.ne.ipl2)   go to 63
          ipti=ipt(itt3-1)
          if(ipti.ne.ipl1.and.ipti.ne.ipl2)   go to 63
          ipti=ipt(itt3)
c checks whether the exchange is necessary.
   63     if(maxmn(xd,yd,ip1,ipti,ipl1,ipl2).eq.0)     go to 64
c modifies the ipt array when necessary.
          ipt(itt3-2)=ipti
          ipt(itt3-1)=ipl1
          ipt(itt3)  =ip1
          ipt(ntt3-1)=ipti
          if(jp2.eq.jpmx)      ipl(jp2t3)=it
          if(jp2.eq.nl0.and.ipl(3).eq.it)     ipl(3)=nt0
c sets flags in the iwl array.
          jwl=jwl+4
          iwl(jwl-3)=ipl1
          iwl(jwl-2)=ipti
          iwl(jwl-1)=ipti
          iwl(jwl)  =ipl2
   64   continue
        nl0=nln
        nlt3=nlnt3
        nlf=jwl/2
        if(nlf.eq.0)      go to 79
c improves the triangulation.
        ntt3p3=ntt3+3
        do 78  irep=1,nrep
          do 76  ilf=1,nlf
            ilft2=ilf*2
            ipl1=iwl(ilft2-1)
            ipl2=iwl(ilft2)
c locates in the ipt array two triangles on both sides of the
c flagged line segment.
            ntf=0
            do 71  itt3r=3,ntt3,3
              itt3=ntt3p3-itt3r
              ipt1=ipt(itt3-2)
              ipt2=ipt(itt3-1)
              ipt3=ipt(itt3)
              if(ipl1.ne.ipt1.and.ipl1.ne.ipt2.and.
     1           ipl1.ne.ipt3)      go to 71
              if(ipl2.ne.ipt1.and.ipl2.ne.ipt2.and.
     1           ipl2.ne.ipt3)      go to 71
              ntf=ntf+1
              itf(ntf)=itt3/3
              if(ntf.eq.2)     go to 72
   71       continue
            if(ntf.lt.2)       go to 76
c determines the vertexes of the triangles that do not lie on the
c line segment.
   72       it1t3=itf(1)*3
            ipti1=ipt(it1t3-2)
            if(ipti1.ne.ipl1.and.ipti1.ne.ipl2)    go to 73
            ipti1=ipt(it1t3-1)
            if(ipti1.ne.ipl1.and.ipti1.ne.ipl2)    go to 73
            ipti1=ipt(it1t3)
   73       it2t3=itf(2)*3
            ipti2=ipt(it2t3-2)
            if(ipti2.ne.ipl1.and.ipti2.ne.ipl2)    go to 74
            ipti2=ipt(it2t3-1)
            if(ipti2.ne.ipl1.and.ipti2.ne.ipl2)    go to 74
            ipti2=ipt(it2t3)
c checks whether the exchange is necessary.
   74 if(maxmn(xd,yd,ipti1,ipti2,ipl1,ipl2).eq.0)
     1         go to 76
c modifies the ipt array when necessary.
            ipt(it1t3-2)=ipti1
            ipt(it1t3-1)=ipti2
            ipt(it1t3)  =ipl1
            ipt(it2t3-2)=ipti2
            ipt(it2t3-1)=ipti1
            ipt(it2t3)  =ipl2
c sets new flags.
            jwl=jwl+8
            iwl(jwl-7)=ipl1
            iwl(jwl-6)=ipti1
            iwl(jwl-5)=ipti1
            iwl(jwl-4)=ipl2
            iwl(jwl-3)=ipl2
            iwl(jwl-2)=ipti2
            iwl(jwl-1)=ipti2
            iwl(jwl)  =ipl1
            do 75  jlt3=3,nlt3,3
              iplj1=ipl(jlt3-2)
              iplj2=ipl(jlt3-1)
              if((iplj1.eq.ipl1.and.iplj2.eq.ipti2).or.
     1           (iplj2.eq.ipl1.and.iplj1.eq.ipti2))
     2                         ipl(jlt3)=itf(1)
              if((iplj1.eq.ipl2.and.iplj2.eq.ipti1).or.
     1           (iplj2.eq.ipl2.and.iplj1.eq.ipti1))
     2                         ipl(jlt3)=itf(2)
   75       continue
   76     continue
          nlfc=nlf
          nlf=jwl/2
          if(nlf.eq.nlfc)      go to 79
c resets the iwl array for the next round.
          jwl=0
          jwl1mn=(nlfc+1)*2
          nlft2=nlf*2
          do 77  jwl1=jwl1mn,nlft2,2
            jwl=jwl+2
            iwl(jwl-1)=iwl(jwl1-1)
            iwl(jwl)  =iwl(jwl1)
   77     continue
          nlf=jwl/2
   78   continue
   79 continue
c rearranges the ipt array so that the vertexes of each triangle are
c listed counter-clockwise.
      do 81  itt3=3,ntt3,3
        ip1=ipt(itt3-2)
        ip2=ipt(itt3-1)
        ip3=ipt(itt3)
        if(side(xd(ip1),yd(ip1),xd(ip2),yd(ip2),xd(ip3),yd(ip3))
     1       .ge.0.0d0)     go to 81
        ipt(itt3-2)=ip2
        ipt(itt3-1)=ip1
   81 continue
      nt=nt0
      return
c
c error exit.
   90 write (nout,2090) ndp0
      go to 93
   91 write (nout,2091) ndp0,ip1,ip2,x1,y1
      go to 93
   92 write (nout,2092) ndp0
   93 write (nout,2093)
      nt=0
      stop
c format statements.
 2090 format(1x/25h ***   ndp less than 4.d0/7h  ndp =,i5)
 2091 format(1x/29h ***   identical data points./
     1       7h  ndp =,i5,5x,5hip1 =,i5,5x,5hip2 =,
     2       i5,5x,4hxd =,e12.4,5x,4hyd =,e12.4)
 2092 format(1x/33h ***   all collinear data points./
     1       7h  ndp =,i5)
 2093 format(35h error detected in routine    ctang)
      end
