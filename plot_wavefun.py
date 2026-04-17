import os
import sys
import numpy as np
import matplotlib.pyplot as plt


# =========================
# Loader
# =========================
def load_complex_data(filepath):
    data = np.loadtxt(filepath)
    x = data[:, 0]
    re = data[:, 1]
    im = data[:, 2]
    return x, re, im


# =========================
# Derived quantities
# =========================
def compute_amplitude(re, im):
    return np.sqrt(re**2 + im**2)


def compute_phase(re, im):
    phase = np.arctan2(im, re)
    return np.unwrap(phase)


# =========================
# Plotting
# =========================
def plot_wavefunction(x, re, im, mode="all", title=""):

    amp = compute_amplitude(re, im)

    if mode == "simple":
        # --- only Re / Im
        plt.figure(figsize=(8, 5))
        plt.plot(x, re, label="Re(ψ)")
        plt.plot(x, im, label="Im(ψ)")
        plt.xlabel(r"$x$")
        plt.ylabel(r"$\psi(x)$")
        plt.legend()
        plt.grid()

    elif mode == "all":
        # --- Re/Im + amplitude
        fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

        axs[0].plot(x, re, label="Re(ψ)")
        axs[0].plot(x, im, label="Im(ψ)")
        axs[0].set_ylabel(r"$\psi(x)$")
        axs[0].legend()
        axs[0].grid()

        axs[1].plot(x, amp, label="|ψ|", linestyle="--")
        axs[1].set_ylabel(r"$|\psi(x)|$")
        axs[1].set_xlabel(r"$x$")
        axs[1].grid()

    elif mode == "phase":
        # --- amplitude + phase
        phase = compute_phase(re, im)

        fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

        axs[0].plot(x, amp)
        axs[0].set_ylabel(r"$|\psi(x)|$")
        axs[0].grid()

        axs[1].plot(x, phase)
        axs[1].set_ylabel(r"$\theta(x)$")
        axs[1].set_xlabel(r"$x$")
        axs[1].grid()

    else:
        raise ValueError("Unknown mode")

    plt.suptitle(title)
    plt.tight_layout()


# =========================
# Main
# =========================
def main():

    # ---- CLI arguments ----
    if len(sys.argv) > 1:
        filepath = sys.argv[1]
    else:
        filepath = "wf_input.dat"

    if len(sys.argv) > 2:
        mode = sys.argv[2]
    else:
        mode = "all"   # options: simple / all / phase

    # ---- load ----
    if not os.path.exists(filepath):
        print("File not found:", filepath)
        return

    x, re, im = load_complex_data(filepath)

    # ---- plot ----
    title = os.path.basename(filepath)
    plot_wavefunction(x, re, im, mode=mode, title=title)

    # ---- save ----
    output_dir = "figures"
    os.makedirs(output_dir, exist_ok=True)

    name = os.path.splitext(os.path.basename(filepath))[0]
    output_path = os.path.join(output_dir, f"{name}_{mode}.png")

    plt.savefig(output_path, dpi=300)
    plt.close()

    print("Saved:", output_path)


# =========================
# Entry point
# =========================
if __name__ == "__main__":
    main()
