from netCDF4 import Dataset as nc
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


def setstr(str):
    if str == 'l':
        return ['lath', 'lonh', 'degrees']
    elif str == 'x':
        return ['yh', 'xh', 'meter']


def create_grd(grdname, ny, nx):

    dimx = 'xi'
    dimy = 'eta'
    ff = nc(grdname, 'w')
    ff.createDimension(dimy, ny)
    ff.createDimension(dimx, nx)
    ff.createDimension('scalar', 1)

    ff.createVariable('h', 'd', (dimy, dimx, ))

    ff.createVariable('xl', 'd', ('scalar',))
    ff.createVariable('el', 'd', ('scalar',))
    ff.createVariable('spherical', 'i', ('scalar',))
    ff.createVariable('f', 'd', (dimy, dimx, ))
    ff.createVariable('pm', 'd', (dimy, dimx, ))
    ff.createVariable('pn', 'd', (dimy, dimx, ))
    ff.createVariable('dmde', 'd', (dimy, dimx, ))
    ff.createVariable('dndx', 'd', (dimy, dimx, ))
    ff.createVariable('angle', 'd', (dimy, dimx, ))
    ff.createVariable('mask_rho', 'd', (dimy, dimx, ))

    ff.close()


def create_ini(grdname, nz, ny, nx):

    # dimt = 'ocean_time'
    dimx = 'xi'
    dimy = 'eta'
    dimz = 'k'
    ff = nc(grdname, 'w')
    ff.createDimension(dimz, nz)
    ff.createDimension(dimy, ny)
    ff.createDimension(dimx, nx)
    ff.createDimension('scalar', 1)
    ff.createVariable('ocean_time', 'd', ('scalar',))

    ff.createVariable('zeta', 'd', (dimy, dimx, ))
    ff.createVariable('ubar', 'd', (dimy, dimx, ))
    ff.createVariable('vbar', 'd', (dimy, dimx, ))
    ff.createVariable('u', 'd', (dimz, dimy, dimx, ))
    ff.createVariable('v', 'd', (dimz, dimy, dimx, ))
    ff.createVariable('temp', 'd', (dimz, dimy, dimx, ))
    ff.createVariable('salt', 'd', (dimz, dimy, dimx, ))

    ff.close()


def create_bry(grdname, nt, nz, ny, nx):

    # dimt = 'ocean_time'
    dimx = 'xi'
    dimy = 'eta'
    dimz = 'k'
    ff = nc(grdname, 'w')
    ff.createDimension('zeta_time', nt)
    ff.createDimension('v2d_time', nt)
    ff.createDimension('v3d_time', nt)
    ff.createDimension('salt_time', nt)
    ff.createDimension('temp_time', nt)
    ff.createDimension(dimz, nz)
    ff.createDimension(dimy, ny)
    ff.createDimension(dimx, nx)
    ff.createVariable('zeta_time', 'd', ('zeta_time', ))
    ff.createVariable('zeta_west', 'd', ('zeta_time', dimy, ))
    ff.createVariable('zeta_east', 'd', ('zeta_time', dimy, ))
    ff.createVariable('zeta_south', 'd', ('zeta_time', dimx, ))
    ff.createVariable('zeta_north', 'd', ('zeta_time', dimx, ))

    ff.createVariable('v2d_time', 'd', ('v2d_time', ))
    ff.createVariable('ubar_west', 'd', ('v2d_time', dimy, ))
    ff.createVariable('ubar_east', 'd', ('v2d_time', dimy, ))
    ff.createVariable('ubar_south', 'd', ('v2d_time', dimx, ))
    ff.createVariable('ubar_north', 'd', ('v2d_time', dimx, ))
    ff.createVariable('vbar_west', 'd', ('v2d_time', dimy, ))
    ff.createVariable('vbar_east', 'd', ('v2d_time', dimy, ))
    ff.createVariable('vbar_south', 'd', ('v2d_time', dimx, ))
    ff.createVariable('vbar_north', 'd', ('v2d_time', dimx, ))

    ff.createVariable('v3d_time', 'd', ('v3d_time', ))
    ff.createVariable('u_west', 'd', ('v3d_time', dimz, dimy, ))
    ff.createVariable('u_east', 'd', ('v3d_time', dimz, dimy, ))
    ff.createVariable('u_south', 'd', ('v3d_time', dimz, dimx, ))
    ff.createVariable('u_north', 'd', ('v3d_time', dimz, dimx, ))
    ff.createVariable('v_west', 'd', ('v3d_time', dimz, dimy, ))
    ff.createVariable('v_east', 'd', ('v3d_time', dimz, dimy, ))
    ff.createVariable('v_south', 'd', ('v3d_time', dimz, dimx, ))
    ff.createVariable('v_north', 'd', ('v3d_time', dimz, dimx, ))

    ff.createVariable('temp_time', 'd', ('temp_time', ))
    ff.createVariable('temp_west', 'd', ('temp_time', dimz, dimy, ))
    ff.createVariable('temp_east', 'd', ('temp_time', dimz, dimy, ))
    ff.createVariable('temp_south', 'd', ('temp_time', dimz, dimx, ))
    ff.createVariable('temp_north', 'd', ('temp_time', dimz, dimx, ))

    ff.createVariable('salt_time', 'd', ('salt_time', ))
    ff.createVariable('salt_west', 'd', ('salt_time', dimz, dimy, ))
    ff.createVariable('salt_east', 'd', ('salt_time', dimz, dimy, ))
    ff.createVariable('salt_south', 'd', ('salt_time', dimz, dimx, ))
    ff.createVariable('salt_north', 'd', ('salt_time', dimz, dimx, ))
    # ff.createVariable('ubar', 'd', (dimt, dimy, dimx, ))
    # ff.createVariable('vbar', 'd', (dimt, dimy, dimx, ))
    # ff.createVariable('u', 'd', (dimt, dimz, dimy, dimx, ))
    # ff.createVariable('v', 'd', (dimt, dimz, dimy, dimx, ))
    # ff.createVariable('temp', 'd', (dimt, dimz, dimy, dimx, ))
    # ff.createVariable('salt', 'd', (dimt, dimz, dimy, dimx, ))

    ff.close()


def create_frc(grdname, nt, ny, nx):

    ff = nc(grdname, 'w')

    ff.createDimension('sms_time', nt)
    ff.createDimension('shf_time', nt)
    ff.createDimension('swf_time', nt)
    ff.createDimension('xi', nx)
    ff.createDimension('xi_u', nx - 1)
    ff.createDimension('eta', ny)
    ff.createDimension('eta_v', ny - 1)

    ff.createVariable('sms_time', 'd', ('sms_time', ))
    ff.createVariable('shf_time', 'd', ('shf_time', ))
    ff.createVariable('swf_time', 'd', ('swf_time', ))

    ff.createVariable('sustr', 'd', ('sms_time', 'eta', 'xi_u', ))
    ff.createVariable('svstr', 'd', ('sms_time', 'eta_v', 'xi', ))
    ff.createVariable('shflux', 'd', ('shf_time', 'eta', 'xi', ))
    ff.createVariable('swflux', 'd', ('swf_time', 'eta', 'xi', ))

    ff.close()


def create_init(initname, grdname, str):
    dimy, dimx, units = setstr(str)

    grid = nc(grdname, 'r')
    lon = grid.variables[dimx][:]
    lat = grid.variables[dimy][:]
    grid.close()

    ff = nc(initname, 'w', format='NETCDF3_CLASSIC')

    ff.createDimension(dimy, np.size(lat))
    ff.createDimension(dimx, np.size(lon))
    ff.createDimension('interface', 3)
    ff.createDimension('LAYER', 2)

    ff.createVariable(dimy, 'd', (dimy, ))
    ff.createVariable(dimx, 'd', (dimx, ))
    ff.createVariable('interface', 'd', ('interface', ))
    ff.createVariable('LAYER', 'd', ('LAYER', ))
    ff.createVariable('ETA', 'd', ('interface', dimy, dimx, ))
    ff.createVariable('u', 'd', ('LAYER', dimy, dimx, ))
    ff.createVariable('v', 'd', ('LAYER', dimy, dimx, ))

    ff.close()
    # copy lon/lat from grid
    ff = nc(initname, 'a')
    ff.variables[dimx][:] = lon
    ff.variables[dimy][:] = lat
    ff.close()


def create_bdy(bdyname, grdname, str):
    dimy, dimx, units = setstr(str)

    grid = nc(grdname, 'r')
    lon = grid.variables[dimx][:]
    lat = grid.variables[dimy][:]
    grid.close()

    ff = nc(bdyname, 'w', format='NETCDF3_CLASSIC')

    ff.createDimension(dimy, np.size(lat))
    ff.createDimension(dimx, np.size(lon))
    ff.createDimension('interface', 3)
    ff.createDimension('LAYER', 2)

    ff.createVariable(dimy, 'd', (dimy, ))
    ff.createVariable(dimx, 'd', (dimx, ))
    ff.createVariable('interface', 'd', ('interface', ))
    ff.createVariable('LAYER', 'd', ('LAYER', ))
    ff.createVariable('ETA', 'd', ('interface', dimy, dimx, ))
    ff.createVariable('u', 'd', ('LAYER', dimy, dimx, ))
    ff.createVariable('v', 'd', ('LAYER', dimy, dimx, ))

    ff.close()

    # copy lon/lat from grid
    ff = nc(bdyname, 'a')
    ff.variables[dimx][:] = lon
    ff.variables[dimy][:] = lat
    ff.close()
#


def create_spng(spngname, grdname, str):
    dimy, dimx, units = setstr(str)

    grid = nc(grdname, 'r')
    lon = grid.variables[dimx][:]
    lat = grid.variables[dimy][:]
    grid.close()

    ff = nc(spngname, 'w', format='NETCDF3_CLASSIC')

    ff.createDimension(dimy, np.size(lat))
    ff.createDimension(dimx, np.size(lon))

    ff.createVariable(dimy, 'd', (dimy, ))
    ff.createVariable(dimx, 'd', (dimx, ))
    ff.createVariable('Idamp', 'd', (dimy, dimx, ))

    ff.close()

    ff = nc(spngname, 'a')
    ff.variables[dimx][:] = lon
    ff.variables[dimy][:] = lat
    ff.close()


######################################################

def hanning_smoother(h):
    [M, L] = np.array(h.shape)
    L = L - 1
    M = M - 1
    Mm = M - 1
    Mmm = M - 2
    Lm = L - 1
    Lmm = L - 2

    h[1:Mm, 1:Lm] = 0.125 * (h[0:Mmm, 1:Lm] + h[2:M, 1:Lm] +
                             h[1:Mm, 0:Lmm] + h[1:Mm, 2:L] +
                             4 * h[1:Mm, 1:Lm])

    h[0, :] = h[1, :]
    h[M, :] = h[Mm, :]
    h[:, 0] = h[:, 1]
    h[:, L] = h[:, Lm]
    return h
