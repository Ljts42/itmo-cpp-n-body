#include "nbody.h"

#include <iostream>

int main()
{
    std::string filename = "../test/etc/planets.txt";
    FastPositionTracker universe(filename);
    std::vector<std::string> planets = {"Earth", "Mars", "Mercury", "Sun", "Venus"};
    for (const auto & planet : planets) {
        Track track = universe.track(planet, 100, 1);
        std::cout << planet << '\t' << track.back().x << ' ' << track.back().y << '\n';
    }
}
