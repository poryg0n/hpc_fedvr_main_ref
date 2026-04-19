import sys

def run_plot(prefix="density", scale="linear", keyword=None):
    if keyword is not None:
        sys.argv = ["script", prefix, scale, keyword]
    else:
        sys.argv = ["script", prefix, scale]

    exec(open('gpt_assisted_plotter.py').read())


if __name__ == "__main__":
    prefix  = sys.argv[1] if len(sys.argv) > 1 else "density"
    scale   = sys.argv[2] if len(sys.argv) > 2 else "linear"
    keyword = sys.argv[3] if len(sys.argv) > 3 else None

    run_plot(prefix, scale, keyword)
