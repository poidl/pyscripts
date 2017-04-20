from netCDF4 import Dataset as nc
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


def setstr(str):
    if str == 'l':
        return ['lath', 'lonh', 'degrees']
    elif str == 'x':
        return ['yh', 'xh', 'meter']


# def create_grd(grdname, lat, lon, str):


def create_grd(grdname, ny, nx):
    # dimy, dimx, units = setstr(str)
    # ff = nc(grdname, 'w', format='NETCDF3_CLASSIC')
    # ff.createDimension(dimy, np.size(lat))
    # ff.createDimension(dimx, np.size(lon))
    dimx = 'xi'
    dimy = 'eta'
    ff = nc(grdname, 'w')
    ff.createDimension(dimy, ny)
    ff.createDimension(dimx, nx)
    ff.createDimension('scalar', 1)

    # ff.createVariable(dimy, 'd', (dimy, ))
    # tmp = ff.variables[dimy]
    # if str == 'l':
    #     setattr(tmp, "long_name", "Latitude")
    # elif str == 'x':
    #     setattr(tmp, "long_name", "y dist")
    # setattr(tmp, "units", units)

    # ff.createVariable(dimx, 'd', (dimx, ))
    # tmp = ff.variables[dimx]
    # if str == 'l':
    #     setattr(tmp, "long_name", "Longitude")
    # elif str == 'x':
    #     setattr(tmp, "long_name", "x dist")
    # setattr(tmp, "units", units)

    # ff.createVariable('D', 'd', (dimy, dimx, ))
    # tmp = ff.variables['D']
    # setattr(tmp, "long_name", "Depth")
    # setattr(tmp, "units", "meter")

    ff.createVariable('h', 'd', (dimy, dimx, ))
    # tmp = ff.variables['D']
    # setattr(tmp, "long_name", "Depth")
    # setattr(tmp, "units", "meter")
    ff.createVariable('xl', 'd', ('scalar',))
    ff.createVariable('el', 'd', ('scalar',))
    ff.createVariable('spherical', 'i', ('scalar',))

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
