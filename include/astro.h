#ifndef ASTRO_H
#define ASTRO_H

#define UA 149597870691.0
#define MASS_SUN 1.98847e30
#define DAY2SEC 86400.0
#define MU_SUN 1.32712440018e+20
#define CORRECTION 0.827363

#define MERCURY 1
#define VENUS 2
#define EARTH 3
#define MARS 4
#define JUPITER 5
#define SATURN 6
#define URANUS 7
#define NEPTUNE 8

#define PI 3.14159265359
#define EXP 2.71828182846

namespace kepler_orbits
{
    // JPL approximate ephemeris of plantes 1800 AD to 2050 AD
    static const float mercury[6] = {0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593};
    static const float venus[6] = {0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255};
    static const float earth[6] = {1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0};
    static const float mars[6] = {1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891};
    static const float jupiter[6] = {5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909};
    static const float saturn[6] = {9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448};
    static const float uranus[6] = {19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630, 74.01692503};
    static const float neptune[6] = {30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574};

    // Cy values
    static const float mercury_cy[6] = {0.00000037, 0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081};
    static const float venus_cy[6] = {0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418};
    static const float earth_cy[6] = {0.00000562, -0.00004392,  -0.01294668, 35999.37244981, 0.32327364, 0.0};
    static const float mars_cy[6] = {0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343};
    static const float jupiter_cy[6] = {-0.00011607, -0.00013253, -0.00183714, 3034.74612775, 0.21252668, 0.20469106};
    static const float saturn_cy[6] = {-0.00125060, -0.00050991,  0.00193609, 1222.49362201, -0.41897216, -0.28867794};
    static const float uranus_cy[6] = {-0.00196176, -0.00004397, -0.00242939, 428.48202785, 0.40805281, 0.04240589};
    static const float neptune_cy[6] = {0.00026291, 0.00005105, 0.00035372, 218.45945325, -0.32241464, -0.00508664};
} // namespace kepler_orbits

namespace planet_props{

    struct Properties{
        double mu;
        double rad;
        double mass;
        double sun_dist;
    };

    static const Properties mercury = {2.2032e13, 2439.5e3, 0.33e24, 57.9e9};
    static const Properties venus = {3.24859e14, 6052e3, 4.87e24, 108.2e9};
    static const Properties earth = {3.986004418e14, 6378e3, 5.97e24, 149.6e9};
    static const Properties mars = {4.282837e13, 3396e3, 0.642e24, 228.0e9};
    static const Properties jupiter = {1.26686534e17, 71492e3, 1898e24, 778.5e9};
    static const Properties saturn = {3.7931187e16, 60268e3, 568e24, 1432.0e9};
    static const Properties uranus = {6.836529e15, 25559e3, 86.8e24, 2867.0e9};
    static const Properties neptune = {5.793939e15, 24764e3, 102e24, 4515.0e9};
} // namespace planet_props

struct orbitalParameters{
    float a0, acy; // semi-major axis
    float e0, ecy; // eccentricity
    float I0, Icy; // inlcinations
    float L0, Lcy; // mean longitude
    float long_peri0, long_pericy; // longitude of preiphasis
    float long_node0, long_nodecy; // acending node lingitude
};
    
#endif //ASTRO_H