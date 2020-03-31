#!/usr/bin/env python

"""
DMFT Self-energy and Green's Function Plotter.

This script plots Greens functions and self energies from DMFT calculations.

"""

import argparse
import os

import matplotlib.pyplot as plt


def plotDMFT(args):
    """This function plots the real and imaginary
    parts of
    1. Green's function
    2. Self energy
    3. Analytically continued Self energy
    """

    gf_file = "G_loc.out." + args.iteration
    se_file = "sig.inp." + args.iteration

    # creating directory for ac
    if not os.path.exists("plots"):
        os.makedirs("plots")

    # plotting Green's functions
    with open(gf_file, "r") as f:
        f.readline()
        lines = f.readlines()
        x = [float(line.split()[0]) for line in lines]
        y1_real = [float(line.split()[1]) for line in lines]  # eg
        y1_imag = [float(line.split()[2]) for line in lines]
        y2_real = [float(line.split()[3]) for line in lines]  # t2g
        y2_imag = [float(line.split()[4]) for line in lines]

    plt.figure(1)
    plt.plot(x, y1_real, label="eg")
    plt.plot(x, y2_real, label="t2g")
    plt.title("Real Green function")
    plt.xlabel("$i{\omega_n}$")
    plt.ylabel("Re $G(i{\omega_n})$")
    plt.legend()
    plt.savefig("./plots/Gf_real.png")
    if args.show:
        plt.show()

    plt.figure(2)
    plt.plot(x, y1_imag, label="eg")
    plt.plot(x, y2_imag, label="t2g")
    plt.title("Imaginary Green function")
    plt.legend()
    plt.xlabel("$i{\omega_n}$")
    plt.ylabel("Im $G(i{\omega_n})$")
    plt.savefig("./plots/Gf_imag.png")
    if args.show:
        plt.show()
    f.close()

    # plotting Self-energies
    with open(se_file, "r") as f:
        for i in range(1, 6):
            f.readline()
        lines = f.readlines()
        x = [float(line.split()[0]) for line in lines]
        y1_real = [float(line.split()[1]) for line in lines]  # eg
        y1_imag = [float(line.split()[2]) for line in lines]
        y2_real = [float(line.split()[3]) for line in lines]  # t2g
        y2_imag = [float(line.split()[4]) for line in lines]

    plt.figure(3)
    plt.plot(x, y1_real, label="eg")
    plt.plot(x, y2_real, label="t2g")
    plt.title("Real Self-Energy")
    plt.xlabel("$i{\omega_n}$")
    plt.ylabel("Re $\Sigma(i{\omega_n})$")
    plt.legend()
    plt.savefig("./plots/Selfenergy_real.png")
    if args.show:
        plt.show()

    plt.figure(4)
    plt.plot(x, y1_imag, label="eg")
    plt.plot(x, y2_imag, label="t2g")
    plt.title("Imaginary Self-Energy")
    plt.legend()
    plt.xlabel("$i{\omega_n}$")
    plt.ylabel("Im $\Sigma(i{\omega_n})$")
    plt.savefig("./plots/Selfenergy_imaginary.png")
    if args.show:
        plt.show()
    f.close()

    # plotting analytically continued Self-energies
    with open("./ac/Sig.out", "r") as f:
        f.readline()
        f.readline()
        lines = f.readlines()
        x = [float(line.split()[0]) for line in lines]
        y1_real = [float(line.split()[1]) for line in lines]  # eg
        y1_imag = [float(line.split()[2]) for line in lines]
        y2_real = [float(line.split()[3]) for line in lines]  # t2g
        y2_imag = [float(line.split()[4]) for line in lines]

    plt.figure(5)
    plt.plot(x, y1_real, label="eg")
    plt.plot(x, y2_real, label="t2g")
    plt.title("Real Self-Energy Analytically continued")
    plt.legend()
    # plt.ylim(-0.1,0.1)
    plt.xlabel("$\omega$")
    plt.ylabel("Re $\Sigma(\omega})$")
    plt.savefig("./plots/Self-EnergyAnalyticallycontinued_real.png")
    if args.show:
        plt.show()

    plt.figure(6)
    plt.plot(x, y1_imag, label="eg")
    plt.plot(x, y2_imag, label="t2g")
    plt.title("Imaginary Self-Energy Analytically continued")
    plt.legend()
    # plt.ylim(-0.1,0.1)
    plt.xlabel("$\omega$")
    plt.ylabel("Im $\Sigma(\omega})$")
    plt.savefig("./plots/Self-EnergyAnalyticallycontinued_imaginary.png")
    if args.show:
        plt.show()
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script plots Greens functions and self energies from DMFT calculations."
    )
    parser.add_argument(
        "iteration", type=str, help="Iteration number of G_loc or sig.inp file.",
    )
    parser.add_argument("-show", help="Show plot on screen", action="store_true")
    args = parser.parse_args()
    plotDMFT(args)
