#!/usr/bin/env python3

from collections import namedtuple
import re

import matplotlib.pyplot as plt
import numpy as np


ANG2AU = 0.5292
AU2EV = 27.2114
NDResult = namedtuple(
            "NDResult",
            "xs numerator h adiabats diabats Us thetas",
)


def nddiabatize(couplings, energies, atom_weights, state_weights, epsilon, fthresh,
                step, xs):
    adiabats = energies - energies.min()
    adiabats *= AU2EV

    points = len(couplings)
    # Numerator
    weighted_couplings = atom_weights[:,None] * couplings
    mags = np.linalg.norm(weighted_couplings.reshape(points, -1), axis=1)
    numerator = np.where(mags <= fthresh, np.zeros_like(mags), mags-fthresh)

    # Denominator
    energy_diffs = np.diff(energies, axis=1).flatten()
    h = np.gradient(energy_diffs, step/ANG2AU)
    denominator = epsilon * state_weights * h

    # theta_tans = numerator / denominator
    twothetas = np.arctan2(numerator, denominator)
    thetas = twothetas / 2
    # def get_P(thetas):
        # return np.array(( np.cos(thetas), np.sin(thetas),
                         # -np.sin(thetas), np.cos(thetas))).T.reshape(points, 2, 2)
    # P = get_P(thetas)
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


def plot(ndres):
    xs = ndres.xs

    fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=4, sharex=True)
    ax0.plot(xs, ndres.numerator)
    ax0.set_ylim(0, 1)
    ax0.set_ylabel("$N_12$")

    ax1.plot(xs, ndres.h)
    ax1.set_ylabel("h")
    ax1.set_ylim(-0.075, 0.015)
    ax1.axhline(0, linestyle="--", color="k")

    ax2.plot(xs, ndres.adiabats)
    ax2.set_xlim(1, 8)
    ax2.set_ylim(0, 15)
    ax2.set_ylabel("$\Delta$E / eV")
    ax2.set_title("Adiabats")

    ax3.plot(xs, ndres.diabats, "o-")
    ax3.set_xlim(1, 8)
    ax3.set_ylim(0, 15)
    ax3.set_ylabel("$\Delta$E / eV")
    ax3.set_title("Diabats")
    ax3.set_xlabel("Li-F / Å")

    fig.tight_layout()
    plt.show()
    # fig.savefig("diabatization.pdf")


if __name__ == "__main__":
    run()
