# integration over a volume
# example computing the core mass

# coding: utf-8

from ester import *
from numpy import *
import argparse
import os

def main():
      a=star2d(ARGS.model)

      nc=a.conv
      print('nb of domains in the core',nc)
      
      # x is the integrated parameter
      print("integrating %s (-x to config) over the core" % ARGS.x)
      is_normalized = True
      x = getattr(a, ARGS.x)
      normalized_to = ARGS.x + "c"
      if is_normalized:
            normalized_value = getattr(a, normalized_to)
      else:
            normalized_value = 1
      # df is the quantity to integrate over the volume
      # dv is the infinitesimal volume
      dv = a.r**2 * a.rz
      df = x * dv

      ncore=sum(a.npts[0:nc])

      # Integrateur radial: I(1,nr)
      b=a.I.ravel() # 2D array with 1 column ==> 1D array
      Icheb=b[:ncore]

      df_l=[]
      dv_l=[]
      # We first integrate over theta on every radial grid point
      for i in range(ncore):
            # dot compute the scalar product of these two:
            res = dot(a.It,df[i,:])
            df_l.append(res)
            # we compute along the integrated core volume
            dv_l.append(dot(a.It,dv[i,:]))

      # theta is normalized on the pi/2 and span over only a quarter of the 
      # sphere so we multiply by pi/2 * 4 = 2*pi
      df_a = 2 * pi * array(df_l)
      dv_a = 2 * pi * array(dv_l)

      # Then we integrate over 'r' on the first ncore grid points
      # X is the integrated x over the whole volume
      X = dot(df_a, Icheb)
      Vcore = dot(dv_a, Icheb)

      print("normalized (on r and %s):" % normalized_to, X)
      # a.r is normalized on a.R so to get the real value:
      X = X * a.R**3
      Vcore = Vcore * a.R**3

      # a.rho is normalized on a.rhoc
      X = X * normalized_value
      print("integrated '%s' on the core:" % ARGS.x, X)
      print("core volume:", Vcore)

      X_average = (X / Vcore)
      print("'%s' average on core = %6.5f" % (ARGS.x, X_average))

      print("ratio of <%s>/%s :" % (ARGS.x, normalized_to), (X_average / normalized_value))

def get_args():
      parser = argparse.ArgumentParser()
      parser.add_argument(
            "model",
            type=str,
            help="the path to a 2D ESTER model",
      )
      parser.add_argument(
            "-x",
            type=str,
            default="rho",
            help="the parameter to integrate over the core, assumed normalized to <x>c in ESTER",
      )
      args = parser.parse_args()
      if not os.path.isfile(args.model):
            print("first argument 'model' must be a path to a valid file")
            exit(1)

      return args

if __name__ == '__main__':
      ARGS = get_args()
      main()
