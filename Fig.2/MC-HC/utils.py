#!/usr/bin/python3


def heatmap(file, limit, scale, mask, xsize, ysize, zsize, normalize):
    import numpy as np
    import pandas as pd

    axis = np.arange(-limit, limit + 1)
    x, y, z = np.meshgrid(axis, axis, axis)
    ball = np.heaviside((limit ** 2 - (x - limit) ** 2 - (y - limit) ** 2 - (z - limit) ** 2), 1)
    output = np.zeros((zsize, ysize, xsize))
    temp = pd.read_csv(file)
    for index, row in temp.iterrows():
        cx = np.rint(row['x'] * scale).astype(int)
        cy = np.rint(row['y'] * scale).astype(int)
        cz = np.rint(row['z'] * scale).astype(int)
        mx = cx - limit
        my = cy - limit
        mz = cz - limit
        Mx = cx + limit + 1
        My = cy + limit + 1
        Mz = cz + limit + 1
        if (mx >= 0) & (my >= 0) & (mz >= 0) & (Mx < xsize) & (My < ysize) & (Mz < zsize):
            if mask[cz, cy, cx] == 1:
                output[mz:Mz, my:My, mx:Mx] += ball

    if normalize:
        output = (output * 50000) / temp.shape[0]
    output = output * mask

    return output


def read_heat_folder(folder, xsize, ysize, zsize):
    import numpy as np
    from os import listdir
    from os.path import join
    import tifffile as tiff

    lista = listdir(folder)
    array = np.zeros((len(lista), zsize, ysize, xsize))
    n = 0
    for name in lista:
        file = join(folder, name)
        array[n, ...] = tiff.imread(file)
        n = n+1

    return array


def ttest(mean1, mean2, std1, std2, n1, n2):
    import numpy as np
    from scipy import stats

    t = (mean1 - mean2)/np.sqrt((std1**2)/n1 + (std2**2)/n2)
    df = ((std1**2)/n1 + (std2**2)/n2)**2/(((std1**2)/n1)**2/(n1-1) +((std2**2)/n2)**2/(n2-1))

    return 2 - 2*stats.t.cdf(t, df)


def pls(x, y):
    import numpy as np
    diag = np.diagflat((np.ones((1, y.shape[0])) @ y) ** (-1))
    m = diag @ y.transpose() @ x
    r = m - np.ones((m.shape[0], 1)) @ ((np.ones((1, m.shape[0])) @ m) / m.shape[0])
    u, s, v = np.linalg.svd(r, full_matrices=False)
    return u, s, v


def procrustes(u, s, v, u0):
    import numpy as np
    n, o, p = np.linalg.svd(np.matmul(u0.transpose(), u), full_matrices=False)
    q = n @ p.transpose()
    vr = v.transpose() @ q
    ur = u @ np.diagflat(s) @ q
    return ur, vr.transpose()


#!/usr/bin/python3


def heatmap(file, limit, scale, mask, xsize, ysize, zsize, normalize):
    import numpy as np
    import pandas as pd

    axis = np.arange(-limit, limit + 1)
    x, y, z = np.meshgrid(axis, axis, axis)
    ball = np.heaviside((limit ** 2 - (x - limit) ** 2 - (y - limit) ** 2 - (z - limit) ** 2), 1)
    output = np.zeros((zsize, ysize, xsize))
    temp = pd.read_csv(file)
    for index, row in temp.iterrows():
        cx = np.rint(row['x'] * scale).astype(int)
        cy = np.rint(row['y'] * scale).astype(int)
        cz = np.rint(row['z'] * scale).astype(int)
        mx = cx - limit
        my = cy - limit
        mz = cz - limit
        Mx = cx + limit + 1
        My = cy + limit + 1
        Mz = cz + limit + 1
        if (mx >= 0) & (my >= 0) & (mz >= 0) & (Mx < xsize) & (My < ysize) & (Mz < zsize):
            if mask[cz, cy, cx] == 1:
                output[mz:Mz, my:My, mx:Mx] += ball

    if normalize:
        output = (output * 50000) / temp.shape[0]
    output = output * mask

    return output


def read_heat_folder(folder, xsize, ysize, zsize):
    import numpy as np
    from os import listdir
    from os.path import join
    import tifffile as tiff

    lista = listdir(folder)
    array = np.zeros((len(lista), zsize, ysize, xsize))
    n = 0
    for name in lista:
        file = join(folder, name)
        array[n, ...] = tiff.imread(file)
        n = n+1

    return array


def ttest(mean1, mean2, std1, std2, n1, n2):
    import numpy as np
    from scipy import stats

    t = (mean1 - mean2)/np.sqrt((std1**2)/n1 + (std2**2)/n2)
    df = ((std1**2)/n1 + (std2**2)/n2)**2/(((std1**2)/n1)**2/(n1-1) +((std2**2)/n2)**2/(n2-1))

    return 2 - 2*stats.t.cdf(t, df)


def pls(x, y):
    import numpy as np
    diag = np.diagflat((np.ones((1, y.shape[0])) @ y) ** (-1))
    m = diag @ y.transpose() @ x
    r = m - np.ones((m.shape[0], 1)) @ ((np.ones((1, m.shape[0])) @ m) / m.shape[0])
    u, s, v = np.linalg.svd(r, full_matrices=False)
    return u, s, v


def procrustes(u, s, v, u0):
    import numpy as np
    n, o, p = np.linalg.svd(np.matmul(u0.transpose(), u), full_matrices=False)
    q = n @ p.transpose()
    vr = v.transpose() @ q
    ur = u @ np.diagflat(s) @ q
    return ur, vr.transpose()


def bootstrap_test(x, y, v0, u0, n, proc):
    import numpy as np
    vdist = np.zeros((n,) + v0.shape)
    m = x.shape[0]
    for i in np.arange(n):
        # generate random index sequence for bootstrapping (i.e. sampling with replacement)
        while True:
            idx = np.random.randint(0, m, m)
            # extract resampled arrays
            xsh = x[idx]
            ysh = y[idx]
            if not np.any(np.all(ysh[..., :] == 0, axis=0)):
                break
        u, s, v = pls(xsh, ysh)
        if proc:
            ur, vr = procrustes(u, s, v, u0)
        else:
            vr = v
        vdist[i, ...] = vr
    vs = np.std(vdist, axis=0)
    return vs

