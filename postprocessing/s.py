import matplotlib.pyplot as plt
import numpy as np

try:
    iti=int(raw_input('first it: '))
    itl=int(raw_input('Last it: '))
except ValueError:
    print "Not a number"

f = open('entropie.txt', 'r')

a = [[], [], [], []] # A list with 3 empty lists
for line in f:
      data = line.split()
      if (data[0] == 'entropy') :
         k=int(data[3])
#	 print k
      else:
         for i, value in enumerate(line.split()):
             a[i].append(float(data[i]))
         a[3].append(k)

g = open('R.txt', 'r')
R = [[], [], []]
for line in g:
      data = line.split()
      if (data[0] == 'it=') :
         k=int(data[1])
#	 print k
      else:
         for i, value in enumerate(line.split()):
             R[i].append(float(data[i]))
         R[2].append(k)


plt.close()
ymm=0.
it=range(iti,itl)
for k_choice in it :
 b = [[], []]
 r = []
 for i in range(k*360):
   if (a[3][i]  == k_choice):
         b[0].append(a[1][i])
         b[1].append(a[2][i])
 for i in range(k*13):
   if (R[2][i]  == k_choice):
         r.append(R[1][i])
#	 print R[1][i]


 #ym=min(b[1][250:360])     
 ym=min(b[1][0:360])     
 ymm=max([ymm,max(b[1][:])])
 plt.plot(b[0][:],b[1][:])
 plt.plot([r[9],r[9]],[-100,100],'b-')
 plt.plot([r[10],r[10]],[-100,100],'r-')
 plt.plot([r[11],r[11]],[-100,100],'k-')
 if (k_choice == max(it)-1):
    plt.plot(b[0][:],b[1][:],'bo')
    plt.plot([r[9],r[9]],[-100,100],'b-x',linewidth=3)
    plt.plot([r[10],r[10]],[-100,100],'r-x',linewidth=3)
    plt.plot([r[11],r[11]],[-100,100],'k-x',linewidth=3)
 if (k_choice == max(it)):
     plt.plot(b[0][:],b[1][:],'ro')
     plt.plot([r[9],r[9]],[-100,100],'b--',linewidth=3)
     plt.plot([r[10],r[10]],[-100,100],'r--',linewidth=3)
     plt.plot([r[11],r[11]],[-100,100],'k--',linewidth=3)
 plt.xlim(xmin=0.,xmax=1.)
 plt.ylim(ymin=ym,ymax=ymm)
plt.show()
