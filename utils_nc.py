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
    dimt = 't'
    ff = nc(grdname, 'w')
    ff.createDimension(dimt, nt)
    ff.createDimension(dimz, nz)
    ff.createDimension(dimy, ny)
    ff.createDimension(dimx, nx)
    ocean_time = ff.createVariable('ocean_time', 'd', (dimt, ))
    ff.createVariable('zeta', 'd', (dimt, dimy, dimx, ))
    ff.createVariable('ubar', 'd', (dimt, dimy, dimx, ))
    ff.createVariable('vbar', 'd', (dimt, dimy, dimx, ))
    ff.createVariable('u', 'd', (dimt, dimz, dimy, dimx, ))
    ff.createVariable('v', 'd', (dimt, dimz, dimy, dimx, ))
    ff.createVariable('temp', 'd', (dimt, dimz, dimy, dimx, ))
    ff.createVariable('salt', 'd', (dimt, dimz, dimy, dimx, ))

    ff.close()


def create_frc(grdname, nt, ny, nx):

    dimx = 'xi'
    dimy = 'eta'
    dimt = 'sms_time'
    ff = nc(grdname, 'w')
    ff.createDimension(dimt, nt)
    ff.createDimension(dimy, ny)
    ff.createDimension(dimx, nx)
    ff.createVariable('sms_time', 'd', (dimt, ))
    var = ff.createVariable('sustr', 'd', (dimt, dimy, dimx, ))
    # var.coordinates = 'lon lat'
    ff.createVariable('svstr', 'd', (dimt, dimy, dimx, ))

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
