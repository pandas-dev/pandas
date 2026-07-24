"""
Various transforms used for by the 3D code
"""

import numpy as np

from matplotlib import _api


def world_transformation(xmin, xmax,
                         ymin, ymax,
                         zmin, zmax, pb_aspect=None):
    """
    Produce a matrix that scales homogeneous coords in the specified ranges
    to [0, 1], or [0, pb_aspect[i]] if the plotbox aspect ratio is specified.
    """
    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin
    if pb_aspect is not None:
        ax, ay, az = pb_aspect
        dx /= ax
        dy /= ay
        dz /= az

    return np.array([[1/dx,    0,    0, -xmin/dx],
                     [   0, 1/dy,    0, -ymin/dy],
                     [   0,    0, 1/dz, -zmin/dz],
                     [   0,    0,    0,        1]])


def _rotation_about_vector(v, angle):
    """
    Produce a rotation matrix for an angle in radians about a vector.
    """
    vx, vy, vz = v / np.linalg.norm(v)
    s = np.sin(angle)
    c = np.cos(angle)
    t = 2*np.sin(angle/2)**2  # more numerically stable than t = 1-c

    R = np.array([
        [t*vx*vx + c,    t*vx*vy - vz*s, t*vx*vz + vy*s],
        [t*vy*vx + vz*s, t*vy*vy + c,    t*vy*vz - vx*s],
        [t*vz*vx - vy*s, t*vz*vy + vx*s, t*vz*vz + c]])

    return R


def _view_axes(E, R, V, roll):
    """
    Get the unit viewing axes in data coordinates.

    Parameters
    ----------
    E : 3-element numpy array
        The coordinates of the eye/camera.
    R : 3-element numpy array
        The coordinates of the center of the view box.
    V : 3-element numpy array
        Unit vector in the direction of the vertical axis.
    roll : float
        The roll angle in radians.

    Returns
    -------
    u : 3-element numpy array
        Unit vector pointing towards the right of the screen.
    v : 3-element numpy array
        Unit vector pointing towards the top of the screen.
    w : 3-element numpy array
        Unit vector pointing out of the screen.
    """
    w = (E - R)
    w = w/np.linalg.norm(w)
    u = np.cross(V, w)
    u = u/np.linalg.norm(u)
    v = np.cross(w, u)  # Will be a unit vector

    # Save some computation for the default roll=0
    if roll != 0:
        # A positive rotation of the camera is a negative rotation of the world
        Rroll = _rotation_about_vector(w, -roll)
        u = np.dot(Rroll, u)
        v = np.dot(Rroll, v)
    return u, v, w


def _view_transformation_uvw(u, v, w, E):
    """
    Return the view transformation matrix.

    Parameters
    ----------
    u : 3-element numpy array
        Unit vector pointing towards the right of the screen.
    v : 3-element numpy array
        Unit vector pointing towards the top of the screen.
    w : 3-element numpy array
        Unit vector pointing out of the screen.
    E : 3-element numpy array
        The coordinates of the eye/camera.
    """
    Mr = np.eye(4)
    Mt = np.eye(4)
    Mr[:3, :3] = [u, v, w]
    Mt[:3, -1] = -E
    M = np.dot(Mr, Mt)
    return M


def _persp_transformation(zfront, zback, focal_length):
    e = focal_length
    a = 1  # aspect ratio
    b = (zfront+zback)/(zfront-zback)
    c = -2*(zfront*zback)/(zfront-zback)
    proj_matrix = np.array([[e,   0,  0, 0],
                            [0, e/a,  0, 0],
                            [0,   0,  b, c],
                            [0,   0, -1, 0]])
    return proj_matrix


def _ortho_transformation(zfront, zback):
    # note: w component in the resulting vector will be (zback-zfront), not 1
    a = -(zfront + zback)
    b = -(zfront - zback)
    proj_matrix = np.array([[2, 0,  0, 0],
                            [0, 2,  0, 0],
                            [0, 0, -2, 0],
                            [0, 0,  a, b]])
    return proj_matrix


def _apply_scale_transforms(xs, ys, zs, axes):
    """
    Apply axis scale transforms to 3D coordinates.

    Transforms data coordinates to transformed coordinates (applying log,
    symlog, etc.) for 3D projection. Preserves masked arrays.
    """
    def transform_coord(coord, axis):
        coord = np.asanyarray(coord)
        data = np.ma.getdata(coord).ravel()
        return axis.get_transform().transform(data).reshape(coord.shape)

    xs_scaled = transform_coord(xs, axes.xaxis)
    ys_scaled = transform_coord(ys, axes.yaxis)
    zs_scaled = transform_coord(zs, axes.zaxis)

    # Preserve combined mask from any masked input
    masks = [np.ma.getmask(a) for a in [xs, ys, zs]]
    if any(m is not np.ma.nomask for m in masks):
        combined = np.ma.mask_or(np.ma.mask_or(masks[0], masks[1]), masks[2])
        xs_scaled = np.ma.array(xs_scaled, mask=combined)
        ys_scaled = np.ma.array(ys_scaled, mask=combined)
        zs_scaled = np.ma.array(zs_scaled, mask=combined)

    return xs_scaled, ys_scaled, zs_scaled


def _proj_transform_vec(vec, M):
    vecw = np.dot(M, vec.data)
    ts = vecw[0:3]/vecw[3]
    if np.ma.isMA(vec):
        ts = np.ma.array(ts, mask=vec.mask)
    return ts[0], ts[1], ts[2]


def _scale_proj_transform_vectors(vecs, axes):
    """
    Apply scale transforms and project vectors.

    Parameters
    ----------
    vecs : ... x 3 np.ndarray
        Input vectors.
    axes : Axes3D
        The 3D axes (used for scale transforms and projection matrix).
    """
    result_shape = vecs.shape
    xs, ys, zs = _apply_scale_transforms(
        vecs[..., 0], vecs[..., 1], vecs[..., 2], axes)
    vec = _vec_pad_ones(xs.ravel(), ys.ravel(), zs.ravel())
    product = np.dot(axes.M, vec)
    tvecs = product[:3] / product[3]
    return tvecs.T.reshape(result_shape)


def _proj_transform_vec_clip(vec, M, focal_length):
    vecw = np.dot(M, vec.data)
    txs, tys, tzs = vecw[0:3] / vecw[3]
    if np.isinf(focal_length):  # don't clip orthographic projection
        tis = np.ones(txs.shape, dtype=bool)
    else:
        tis = (-1 <= txs) & (txs <= 1) & (-1 <= tys) & (tys <= 1) & (tzs <= 0)
    if np.ma.isMA(vec[0]):
        tis = tis & ~vec[0].mask
    if np.ma.isMA(vec[1]):
        tis = tis & ~vec[1].mask
    if np.ma.isMA(vec[2]):
        tis = tis & ~vec[2].mask

    txs = np.ma.masked_array(txs, ~tis)
    tys = np.ma.masked_array(tys, ~tis)
    tzs = np.ma.masked_array(tzs, ~tis)
    return txs, tys, tzs, tis


def inv_transform(xs, ys, zs, invM):
    """
    Transform the points by the inverse of the projection matrix, *invM*.
    """
    vec = _vec_pad_ones(xs, ys, zs)
    vecr = np.dot(invM, vec)
    if vecr.shape == (4,):
        vecr = vecr.reshape((4, 1))
    for i in range(vecr.shape[1]):
        if vecr[3][i] != 0:
            vecr[:, i] = vecr[:, i] / vecr[3][i]
    return vecr[0], vecr[1], vecr[2]


def _vec_pad_ones(xs, ys, zs):
    if np.ma.isMA(xs) or np.ma.isMA(ys) or np.ma.isMA(zs):
        return np.ma.array([xs, ys, zs, np.ones_like(xs)])
    else:
        return np.array([xs, ys, zs, np.ones_like(xs)])


def proj_transform(xs, ys, zs, M):
    """
    Transform the points by the projection matrix *M*.
    """
    vec = _vec_pad_ones(xs, ys, zs)
    return _proj_transform_vec(vec, M)


@_api.deprecated("3.10")
def proj_transform_clip(xs, ys, zs, M):
    vec = _vec_pad_ones(xs, ys, zs)
    return _proj_transform_vec_clip(vec, M, focal_length=np.inf)


def _scale_proj_transform_clip(xs, ys, zs, axes):
    """
    Apply scale transforms, project, and return clipping result.

    Returns txs, tys, tzs, tis.
    """
    xs, ys, zs = _apply_scale_transforms(xs, ys, zs, axes)
    vec = _vec_pad_ones(xs, ys, zs)
    return _proj_transform_vec_clip(vec, axes.M, axes._focal_length)


def _proj_trans_points(points, M):
    points = np.asanyarray(points)
    xs, ys, zs = points[:, 0], points[:, 1], points[:, 2]
    return proj_transform(xs, ys, zs, M)


def _scale_proj_transform(xs, ys, zs, axes):
    """
    Apply scale transforms and project.

    Combines `_apply_scale_transforms` and `proj_transform` into a single
    call. Returns txs, tys, tzs.
    """
    xs, ys, zs = _apply_scale_transforms(xs, ys, zs, axes)
    return proj_transform(xs, ys, zs, axes.M)
