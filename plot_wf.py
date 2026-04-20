import numpy as np
import matplotlib.pyplot as plt


def load_wavefunction(filename):
    t = None
    omg = None

    x = []
    re = []
    im = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()

            # --- header parsing
            if line.startswith("#"):
                if "t" in line:
                    try:
                        t = float(line.split("=")[1])
                    except:
                        pass
                if "omg" in line:
                    try:
                        omg = float(line.split("=")[1])
                    except:
                        pass
                continue

            # --- data
            parts = line.split()
            if len(parts) == 3:
                xi, ri, ii = map(float, parts)
                x.append(xi)
                re.append(ri)
                im.append(ii)

    x = np.array(x)
    psi = np.array(re) + 1j * np.array(im)

    return x, psi, t, omg


def plot_wavefunction(x, psi, t=None, omg=None, show=True, save=None):
    density = np.abs(psi)**2

    plt.figure(figsize=(8, 5))

    plt.plot(x, density, label=r"$|\psi(x)|^2$")
    plt.plot(x, psi.real, "--", label="Re(ψ)")
    plt.plot(x, psi.imag, "--", label="Im(ψ)")

    title = "Wavefunction"
    if t is not None:
        title += f"  |  t = {t:.4f}"
    if omg is not None:
        title += f"  |  ω = {omg:.4f}"

    plt.title(title)
    plt.xlabel("x")
    plt.ylabel("Amplitude")
    plt.yscale("log")


    plt.legend()
    plt.grid()

    if save:
        plt.savefig(save, dpi=150)

    if show:
        plt.show()

    plt.close()


# --- quick usage ---
if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python plot_wf.py file.dat")
        sys.exit(1)

    filename = sys.argv[1]

    x, psi, t, omg = load_wavefunction(filename)
    plot_wavefunction(x, psi, t, omg, save="init.png")
