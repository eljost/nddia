from collections import namedtuple

import numpy as np
from nddia.constants import ANG2AU, AU2EV


NDResult = namedtuple(
    "NDResult",
    "xs numerator h adiabats diabats Us thetas",
)


def diabatize(
    couplings, energies, atom_weights, state_weights, epsilon, fthresh, step, xs
):
    adiabats = energies - energies.min()
    adiabats *= AU2EV

    points = len(couplings)

    # Numerator
    weighted_couplings = atom_weights[:, None] * couplings
    mags = np.linalg.norm(weighted_couplings.reshape(points, -1), axis=1)
    numerator = np.where(mags <= fthresh, np.zeros_like(mags), mags - fthresh)

    # Denominator
    energy_diffs = np.diff(energies, axis=1).flatten()
    h = np.gradient(energy_diffs, step / ANG2AU)
    denominator = epsilon * state_weights * h

    twothetas = np.arctan2(numerator, denominator)
    thetas = twothetas / 2
    """
    def get_P(thetas):
        return np.array(( np.cos(thetas), np.sin(thetas),
                         -np.sin(thetas), np.cos(thetas))).T.reshape(points, 2, 2)
    P = get_P(thetas)
    """
    Us = list()
    diabats = list()
    for ens_, theta in zip(energies, thetas):
        cos = np.cos(theta)
        sin = np.sin(theta)
        P = np.array(((cos, sin), (-sin, cos)))
        V = np.diag(ens_)
        U = np.linalg.inv(P).dot(V).dot(P)
        Us.append(U)
        diabats.append(np.diag(U))

    Us = np.array(Us)
    diabats = np.array(diabats)
    diabats -= diabats.min()
    diabats *= AU2EV

    ndres = NDResult(
        xs=xs,
        numerator=numerator,
        h=h,
        adiabats=adiabats,
        diabats=diabats,
        Us=Us,
        thetas=thetas,
    )
    return ndres


if __name__ == "__main__":
    run()
