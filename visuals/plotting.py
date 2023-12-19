import matplotlib.pyplot as plt
import numpy as np
import datetime

import pykep
import pykep as pk
from pykep.core import propagate_lagrangian


AU = 149597870691.0
N = 60
MU = 1.32712440018e+20
SEC2DAY = 1.1574074074074073e-05
DAY2SEC = 86400
JD2000 = 2451544.50000

count = 0


class Planet:
    def __init__(self, name: str, at, coordinates=None, color='gray', orbit_color='gray'):
        self.name = name
        self.at = at
        self.coordinates = coordinates
        self.color = color
        self.orbit_color = orbit_color

        if self.name not in ['mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune']:
            raise Exception(f"Unknown solar system planet named {self.name}")

        if isinstance(self.at, str):
            try:
                self.at = datetime.datetime.strptime(self.at, "%Y-%m-%d")
                self.at = self.at.toordinal() + 1721424.5
            except ValueError:
                raise ValueError(f'"at" attribute has str format but the value is not of the date type: "%Y-%m-%d"')
        elif isinstance(self.at, (float, int)):
            if self.at < 1721424.5:
                raise ValueError('"at" attribute must be a valid Julian Date.')
        else:
            raise ValueError('invalid data type from "at" attribute, must be str, int or float.')

    def plot(self, ax, ax2d=None):
        """Plot the planet's orbit and the location if the coordinates have been passed"""
        plnt = pk.planet.jpl_lp(self.name)
        T = plnt.compute_period(pk.epoch(self.at, julian_date_type="jd")) * SEC2DAY

        # Uranus and Neptune have really long periods above 2050 -> Out of range ephemeris
        if self.at + T >= 2469807.5:
            T = 2469807.5 - self.at - 1

        when = np.linspace(0, T, N)

        x = np.array([0.0] * N)
        y = np.array([0.0] * N)
        z = np.array([0.0] * N)

        for i, day in enumerate(when):
            r, v = plnt.eph(pk.epoch(self.at + day, julian_date_type="jd"))
            x[i] = r[0] / AU
            y[i] = r[1] / AU
            z[i] = r[2] / AU

        # 3D plot
        ax.plot(x, y, z, color=self.orbit_color)

        # 2D plot
        if ax2d is not None:
            ax2d.plot(x, y, color=self.orbit_color)

        # Plot also the planet location if needed
        if self.coordinates is not None:
            ax.scatter(self.coordinates[0]/AU, self.coordinates[1]/AU, self.coordinates[2]/AU,
                       marker='o', alpha=0.8, label=self.name, color=self.color)
        if self.coordinates is not None and ax2d is not None:
            ax2d.scatter(self.coordinates[0]/AU, self.coordinates[1]/AU, marker='o', alpha=0.8, label=self.name,
                         color=self.color)

    def plot_animate(self, ax):
        """
        Plot the planet's orbit for the animation.
        Only plots the orbit in a 2D plot. NO location.
        Ask Iker for further info on this function and all the other animate ones.
        """
        plnt = pk.planet.jpl_lp(self.name)
        T = plnt.compute_period(pk.epoch(self.at, julian_date_type="jd")) * SEC2DAY
        if self.at + T >= 2469807.5:  # Uranus and Neptune have really long periods above 2050 -> Out of range ephemeris
            T = 2469807.5 - self.at - 1

        when = np.linspace(0, T, N)

        x = np.array([0.0] * N)
        y = np.array([0.0] * N)
        z = np.array([0.0] * N)

        for i, day in enumerate(when):
            r, v = plnt.eph(pk.epoch(self.at + day, julian_date_type="jd"))
            x[i] = r[0] / AU
            y[i] = r[1] / AU
            z[i] = r[2] / AU

        ax.plot(x, y, color=self.orbit_color)


class Transfer:
    def __init__(self, planet_origin: Planet, planet_dest: Planet, velocity, color='teal'):
        self.planet_origin = planet_origin
        self.planet_dest = planet_dest
        self.velocity = velocity
        self.color = color
        self.v_evolution = list()

        if not isinstance(self.planet_origin, Planet) or not isinstance(self.planet_dest, Planet):
            raise ValueError(f'Planet origin and destination must be a reference to Planet type object')

    def solve_transfer(self):
        """Computes the transfer using the lagragian propagation. Stores the values where it corresponds"""

    def plot(self, ax, ax2d=None, animate=False, planets=None):
        """Plot the transfer between the planet of origin and the next one."""
        x = np.zeros(N, )
        y = np.zeros(N, )
        z = np.zeros(N, )

        r = self.planet_origin.coordinates
        t = (self.planet_dest.at - self.planet_origin.at) * DAY2SEC
        ts = t / (N - 1)

        v = self.velocity
        for i in range(N):
            x[i] = r[0] / AU
            y[i] = r[1] / AU
            z[i] = r[2] / AU
            r, v = propagate_lagrangian(r, v, ts, MU)
            self.v_evolution.append(np.linalg.norm(v))

        ax.plot(x, y, z, color=self.color)

        if ax2d is not None:
            ax2d.plot(x, y, color=self.color)

    def plot_animate(self, ax, planets, n):
        """
        Plot the transfer orbit and planets moving for the animation.
        Only plots the orbit in a 2D plot. NO location.
        Ask Iker for further info on this function and all the other animate ones.
        """
        global count
        x = list()
        y = list()
        z = list()

        nn = int(planets[n + 1].at - planets[n].at)//5

        r = self.planet_origin.coordinates
        t = (self.planet_dest.at - self.planet_origin.at) * DAY2SEC
        ts = t / (nn - 1)

        v = self.velocity
        # T = np.linspace(planets[n].at, planets[n + 1].at, nn)
        for i in range(nn):
            x.append(r[0] / AU)
            y.append(r[1] / AU)

            r, v = propagate_lagrangian(r, v, ts, MU)

            ax.plot(x, y, color=self.color)

            # Plot the spacecraft
            spacecraft = ax.scatter(x[-1], y[-1], color=self.color)

            # Plot the planets
            points = list()
            T = np.linspace(planets[n].at, planets[n+1].at, nn)
            # at = T[N*n+i]
            for p in planets:
                plnt = pk.planet.jpl_lp(p.name)
                r_eph, v_eph = plnt.eph(pk.epoch(T[i], julian_date_type="jd"))
                points.append(ax.scatter(r_eph[0] / AU, r_eph[1] / AU, color='gray'))

            plt.savefig(f'animate/plot{count}.png')
            count += 1
            # Remove the moving points
            spacecraft.remove()

            for pp in points:
                pp.remove()