# import bcfar_fort, bcwall_fort, halo_fort, math
from bin.Field import mean

def transfer_down(self, model, workspace1, workspace2):
    """Sets the geometry values in the halo
        
        Args:
            model: The physics model
            workspace1: The finer Workspace
            workspace2: The coarser Workspace
        """
    # get padding
    pad = self.padding

    # get geometry dictionary
    geom1 = workspace1.get_geometry()
    geom2 = workspace2.get_geometry()

    # variables
    ev = workspace1.get_field('ev', model.className)
    lv = workspace1.get_field('lv', model.className)
    evc = workspace2.get_field('ev', model.className)
    lvc = workspace2.get_field('lv', model.className)

    # dims
    [nx, ny] = workspace1.field_size()

    # coarse mesh dims
    ratio = 2
    nxc = int(nx/ratio)
    nyc = int(ny/ratio)

    # parameters
    kvis = model.params['kvis']

    #     if (kvis.gt.0) then
    # c
    # c     transfer the molecular and turbulent viscosity to the coarse grid
    # c
    if kvis > 0:

        jj        = 1
        for j in range(pad, ny+pad, 2):
            jj        = jj  +1
            ii        = 1
            for i in range(pad, nx+pad, 2):
                ii        = ii  +1
                lvc[ii, jj] = mean(lv[i:i+2, j:j+2])
                evc[ii, jj] = mean(lv[i:i+2, j:j+2])

        # c
        # c     set the boundary values at i=1 and i=ie
        # c
        jj        = 1
        for j in range(pad, ny+pad, ratio):
            jj        = jj  +1

            lvc[1    ,jj] = mean(lv[1   , j:j+ratio])
            lvc[nxc+pad,jj] = mean(lv[nx+pad, j:j+ratio])
            evc[1    ,jj] = mean(ev[1   , j:j+ratio])
            evc[nxc+pad,jj] = mean(ev[nx+pad, j:j+ratio])

        # c
        # c     set the boundary values at j=1 and j=je
        # c
        ii        = 1
        for i in range(pad, nx+pad, ratio):
            ii        = ii  +1

            lvc[ii,1    ] = mean(lv[i:i+ratio, 1   ])
            lvc[ii,nyc+pad] = mean(lv[i:i+ratio, ny+pad])
            lvc[ii,1    ] = mean(lv[i:i+ratio, 1   ])
            lvc[ii,nyc+pad] = mean(lv[i:i+ratio, ny+pad])
