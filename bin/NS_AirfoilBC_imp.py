import bcfar_fort, bcwall_fort, halo_fort

class NS_AirfoilBC_imp():

    def bc_far(self, this, workspace, state, p):
        # get geometry dictionary
        geom = workspace.get_geom()
        
        # dims
        il = 0
        jl = 0
        ie = 0
        je = 0
        itl = 0
        itu = 0
        
        # flo_var
        w = 0
        rlv = 0
        rev = 0
        
        # mesh_var
        x = 0
        xc = 0
        
        # out_var
        cp = 0
        cf = 0
        
        # flo_param
        gamma = 0
        rm = 0
        rho0 = 0
        p0 = 0
        h0 = 0
        c0 = 0
        u0 = 0
        v0 = 0
        ca = 0
        sa = 0
        re = 0
        prn = 0
        prt = 0
        scal = 0
        chord = 0
        xm = 0
        ym = 0
        kvis = 0
        
        # solv_param
        bc = 0
        
        # mg_param
        mode = 0        
        
        bcfar_fort.bcfar(il, jl, ie, je, itl, itu,
                         w, p, rlv, rev,
                         x, xc, 
                         cp, cf,
                         gamma,rm,rho0,p0,h0,c0,u0,v0,ca,sa,re,prn,prt,scal,chord,xm,ym,kvis,
                         bc,
                         mode)


    def bc_wall(self, this, workspace, state, p, vol, rev):
        # get geometry dictionary
        geom = workspace.get_geom()
        
        # dims
        ny = 0
        il = 0
        ie = 0
        ib = 0
        itl = 0
        itu = 0
        
        # flo_var
        w = 0
        p = 0
        
        # mesh_var
        x = 0
        
        # flo_param
        rm = 0
        sa = 0
        kvis = 0
        
        # solv_param
        isym = 0
        
        bcwall_fort.bcwall(ny, il, ie, ib, itl, itu, 
                           w, p,
                           x,
                           rm, sa, kvis,
                           isym)


    def halo(self, this, workspace, state, p, vol):
        # get geometry dictionary
        geom = workspace.get_geom()
        
        # dims
        il = 0
        jl = 0
        ie = 0
        je = 0
        ib = 0
        jb = 0
        itl = 0
        itu = 0
        
        # flo_var
        w = 0
        # also p
        
        # mesh_var
        x = 0
        # also vol
        
        halo_fort.halo(il, jl, ie, je, ib, jb, itl, itu,
             w, p,
             x, vol)

    
    def transfer_down(self, this, workspace1, workspace2, rev1, rlv1, rev2, rlv2):
        # get geometry dictionary
        geom1 = workspace1.get_geom()
        geom2 = workspace2.get_geom()