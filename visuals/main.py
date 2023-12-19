import matplotlib.pyplot as plt
import numpy as np

import pykep as pk
from pykep.orbit_plots import plot_planet, plot_lambert
from pykep import AU, DAY2SEC

import plotting
import json


def load(visual_js):
    """Load the visualize file"""
    planets = list()
    transfers = list()

    if len(visual_js["Planets"]) != len(visual_js["Transfers"]) + 1:
        raise ValueError(f'The number of transfers and planets does NOT add up.')

    for p in visual_js["Planets"]:
        planets.append(plotting.Planet(p["Name"], p["At"], p["Coordinates"], p["Color"])) #TODO: Repeat p["Color"] again to use the same color for the orbit

    for idx, t in enumerate(visual_js["Transfers"]):
        transfers.append(plotting.Transfer(planets[idx], planets[idx + 1], t["Velocity"], t["Color"]))

    for i in range(1, len(planets)):
        if planets[i-1].at >= planets[i].at:
            raise ValueError(f'"At" time values at planets are not ordered. Incoherent')

    times = visual_js["Extra"]

    return planets, transfers, times


def visualize(planets, transfers):
    fig = plt.figure(figsize=(10, 8))
    ax3d = plt.subplot(1, 1, 1, projection='3d')
    ax3d.set_xlabel('X (AU)')
    ax3d.set_ylabel('Y (AU)')
    ax3d.set_zlabel('Z (AU)')

    # TODO: 2D plot, PASS AS ARGUMENT IN THE PLOT CALL! If not, only 3D plot.
    fig = plt.figure(figsize=(6, 6))
    ax2d = plt.subplot(1, 1, 1)
    ax2d.set_xlabel('X (AU)')
    ax2d.set_ylabel('Y (AU)')
    ax2d.scatter(0, 0, s=50, marker='o', color='yellow', edgecolor='darkorange')

    v = list()
    for p in planets:
        p.plot(ax3d, ax2d)

    for t in transfers:
        t.plot(ax3d, ax2d)
        v.extend(t.v_evolution)
    x = np.linspace(planets[0].at, planets[-1].at, len(v))

    # SUN
    ax3d.scatter(0, 0, 0, s=50, marker='o', color='yellow', edgecolor='darkorange')

    # Speed plot
    fig = plt.figure(figsize=(10, 4))
    ax_v = plt.subplot(1, 1, 1)
    ax_v.set_ylabel("Speed (m/s)")
    ax_v.set_xlabel("Time (JD)")
    ax_v.plot(x, v)

    plt.show()


def animate(planets, transfers, times):
    fig = plt.figure(figsize=(6, 6))
    ax2d = plt.subplot(1, 1, 1)
    ax2d.set_xlabel('X (AU)')
    ax2d.set_ylabel('Y (AU)')
    ax2d.scatter(0, 0, s=50, marker='o', color='yellow', edgecolor='darkorange')

    for p in planets:
        p.plot_animate(ax2d)

    for n, t in enumerate(transfers):
        t.plot_animate(ax2d, planets, n)


def main():
    """
    Main call. All the file parsing and plotting.
    If the script is called from outside, from mga, it will use the visualize.json file.
    Otherwise, if called from the current directory, it will use the auxiliar one.
    THis is done to be able to plot and do things from here directly.
    """
    # Chose the file to use depending on from where the script is called.
    try:
        with open('visuals/visualize.json', 'rb') as f:
            visual_js = json.load(f)
    except FileNotFoundError:
        with open('visualize_aux.json', 'rb') as f:
            visual_js = json.load(f)

    # Load everything and visualize.
    planets, transfers, times = load(visual_js)
    visualize(planets, transfers)

    # As Iker for this animate function is needed. Otherwise, ignore it...
    # It doesn't animate like that, it just creates the frame in a folder, later composed to video using ffmpeg.
    # animate(planets, transfers, times)


if __name__ == '__main__':
    main()