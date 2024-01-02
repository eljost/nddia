import matplotlib.pyplot as plt
import numpy as np

from nddia.constants import AU2EV
from nddia.parser import parse_molcas_couplings
from nddia.main import diabatize


def plot(ndres, show=False):
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
    if show:
        plt.show()
    return fig


def test_lif():
    energies = np.loadtxt("07_foreach_sym.energies")
    ens = energies - energies.min()
    ens *= AU2EV

    fn = "07_foreach_sym.out"
    with open(fn) as handle:
        text = handle.read()
    couplings = parse_molcas_couplings(text)

    w = np.ones(couplings.shape[1])
    step = 0.1
    xs = 0.8 + np.arange(len(couplings)) * step
    nd_kwargs = {
        "couplings": couplings,
        "energies": energies,
        "atom_weights": w,
        "state_weights": 1,
        "epsilon": 160 / 2.0,
        "fthresh": 0,
        "step": step,
        "xs": xs,
    }
    ndres = diabatize(**nd_kwargs)
    fig = plot(ndres)
    fig.savefig("lif_nddia.pdf")
