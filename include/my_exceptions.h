#ifndef MY_EXCEPTIONS_H
#define MY_EXCEPTIONS_H

#include <exception>

struct UnknowPLanet : public std::exception {
    const char* what() const throw () {
        return "Unrecognized planet.";
    }
};

struct PlanetUnreferended : public std::exception {
    const char* what() const throw () {
        return "The planet instance redferenced on the current flyby instance has been destroyed. Nullptr.";
    }
};

struct negativeTime : public std::exception {
    const char* what() const throw () {
        return "Time cannot be negative.";
    }
};

struct distinctFlybyPlanet : public std::exception {
    const char* what() const throw () {
        return "Incoming and departure planet on a flyby are not the same. Transfer must be arriving and departing at the same body.";
    }
};

struct timeRange : public std::exception {
    const char* what() const throw () {
        return "Time ranges are not coherent (min >= max)";
    }
};

#endif //MY_EXCEPTIONS_H