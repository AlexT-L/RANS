# fortran module
import bcfar_fort 

def bc_far(self, model, workspace, state):
    # define helper function for getting fields
    def get(varName):
        return workspace.get_field(varName, model.className)

    # get geometry dictionary
    geom = workspace.get_geometry()
    dims = workspace.get_dims()
    
    # dims
    [nx, ny] = workspace.field_size()
    il = nx+1
    jl = ny+1
    ie = nx+2
    je = ny+2
    ib = nx+3
    jb = nx+3
    itl = dims['itl']
    itu = dims['itu']
    
    # flo_var
    w = state
    lv = get('lv')
    ev = get('ev')
    p = get("p")
    
    # mesh_var
    coords = workspace.get_field("x")
    x = coords
    coords = workspace.get_field("xc")
    xc = coords
    
    # out_var
    cp = workspace.get_field("cp", self.className)
    cp = cp
    cf = workspace.get_field("cf", self.className)
    cf = cf
    
    # flo_param
    gamma = model.params['gamma']
    rm = model.params['rm']
    rho0 = model.params['rho0']
    p0 = model.params['p0']
    h0 = model.params['h0']
    c0 = model.params['co']
    u0 = model.params['u0']
    v0 = model.params['v0']
    ca = model.params['ca']
    sa = model.params['sa']
    re = model.params['re']
    prn = model.params['prn']
    prt = model.params['prt']
    scal = geom['scal']
    chord = geom['chord']
    xm = geom['xm']
    ym = geom['ym']
    kvis = model.params['kvis']
    
    # solv_param
    bc = self.bc
    
    # mg_param
    mode = 1
    if workspace.is_finest():
        mode = 0
    
    bcfar_fort.bcfar(il, jl, ie, je, itl+1, itu+1,
                        w, p, lv, ev,
                        x, xc, 
                        cp, cf,
                        gamma,rm,rho0,p0,h0,c0,u0,v0,ca,sa,re,prn,prt,scal,chord,xm,ym,kvis,
                        bc,
                        mode)
