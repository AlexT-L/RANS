# append to path so we can access files
import sys
sys.path.append("../../../")

from Field import mean

def metric(self):
    # parameters
    [nx, ny] = self.divisions

    # get arrays
    x = self.fields['x']
    xc = self.fields['xc']
    vol = self.fields['vol']

    # calculate cell centers
    for i in range(nx):
        for j in range(ny):
            xc[i,j,0] = mean(x[i:i+2, j:j+2, 0])
            xc[i,j,1] = mean(x[i:i+2, j:j+2, 1])

            xrh = x[i+1,j+1,0] - x[i  ,j,0]
            yrh = x[i+1,j+1,1] - x[i  ,j,1]
            xlh = x[i  ,j+1,0] - x[i+1,j,0]
            ylh = x[i  ,j+1,1] - x[i+1,j,1]

            vol[i,j] = 0.5*(xrh*ylh - yrh*xlh)