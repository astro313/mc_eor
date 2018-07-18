

def setup_plot():
    import matplotlib
    # print(matplotlib.matplotlib_fname())

    matplotlib.rcParams.update({'figure.figsize': (8, 5)    # inches
                                , 'font.size': 10      # points
                                , 'legend.fontsize': 8      # points
                                , 'lines.linewidth': 2       # points
                                , 'axes.linewidth': 1       # points
                                , 'text.usetex': True    # Use LaTeX to layout text
                                , 'font.family': "serif"  # Use serifed fonts
                                , 'xtick.major.size': 10     # length, points
                                , 'xtick.major.width': 1     # points
                                , 'xtick.minor.size': 6     # length, points
                                , 'xtick.minor.width': 1     # points
                                , 'ytick.major.size': 10     # length, points
                                , 'ytick.major.width': 1     # points
                                , 'ytick.minor.size': 6     # length, points

                                , 'ytick.minor.width': 1     # points
                                , 'font.serif': ("Times", "Palatino", "Computer Modern Roman", "New Century Schoolbook", "Bookman"), 'font.sans-serif': ("Helvetica", "Avant Garde", "Computer Modern Sans serif"), 'font.monospace': ("Courier", "Computer Modern Typewriter"), 'font.cursive': "Zapf Chancery"
                                })
    import matplotlib.pyplot as plt
    return plt


def col_f(ii, cm=None):
    if cm is None:
        cm = plt.get_cmap('gist_heat')
    return cm(ii)

