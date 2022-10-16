#pragma once

#include <memory>
#include <string>
#include <vector>

const double G = 6.67e-11;
const double Theta = 0.5;

struct Cartesian
{
    double x = 0;
    double y = 0;

    Cartesian() = default;
    Cartesian(double, double);

    Cartesian & operator+=(const Cartesian &);
};

// Quadrant representation, required for Problem 2
class Quadrant
{
    Cartesian m_center;
    double m_length;

public:
    // Create quadrant with center (x, y) and size 'lenth'
    Quadrant(Cartesian, double);

    // Test if point (x, y) is in the quadrant
    bool contains(const Cartesian &) const;
    double length() const;

    // The four methods below construct new Quadrant representing sub-quadrant of the invoking quadrant
    Quadrant nw() const;
    Quadrant ne() const;
    Quadrant sw() const;
    Quadrant se() const;

    friend std::ostream & operator<<(std::ostream &, const Quadrant &);
};

// Single body representation, required for Problem 1 and Problem 2
class Body
{
    std::string m_name;
    double m_weight;
    Cartesian m_coord;
    Cartesian m_velocity;
    Cartesian m_force;

public:
    Body() = default;
    Body(const Body &) = default;
    Body(const std::string &, double, Cartesian, Cartesian);

    std::string getName() const;
    double getWeight() const;
    Cartesian getCoord() const;
    Cartesian getForse() const;
    Cartesian getSpeed() const;

    double distance(const Body &) const;

    // calculate the force-on current body by the 'b' and add the value to accumulated force value
    void add_force(const Body &);
    // reset accumulated force value
    void reset_force();

    // update body's velocity and position
    void update(double);

    friend std::ostream & operator<<(std::ostream &, const Body &);
    friend std::istream & operator>>(std::istream &, Body &);

    // The methods below to be done for Burnes-Hut algorithm only
    // Test if body is in quadrant
    bool in(const Quadrant) const;
    // Create new body representing center-of-mass of the invoking body and 'b'
    Body plus(const Body &) const;

    ~Body() = default;
};

// Burnes-Hut tree representation, required for Problem 2
class BHTreeNode
{
    double m_total_weight = 0;
    Cartesian m_mass_center{0, 0};
    Quadrant m_borders;
    std::unique_ptr<Body> m_body = nullptr;

    std::unique_ptr<BHTreeNode> m_nw = nullptr;
    std::unique_ptr<BHTreeNode> m_ne = nullptr;
    std::unique_ptr<BHTreeNode> m_sw = nullptr;
    std::unique_ptr<BHTreeNode> m_se = nullptr;

public:
    BHTreeNode(const Quadrant &);

    void insert(const Body &);
    // Update net acting force-on 'b'
    void update_force(Body &);
};

using Track = std::vector<Cartesian>;

class PositionTracker
{
protected:
    std::vector<Body> m_bodies;
    double m_size;
    PositionTracker(const std::string &);

public:
    virtual Track track(const std::string &, std::size_t, std::size_t) = 0;
    virtual ~PositionTracker() = default;
};

class BasicPositionTracker : public PositionTracker
{
public:
    BasicPositionTracker(const std::string &);
    ~BasicPositionTracker() = default;
    Track track(const std::string &, std::size_t, std::size_t) override;
};

class FastPositionTracker : public PositionTracker
{
    BHTreeNode * m_root;

public:
    FastPositionTracker(const std::string &);
    Track track(const std::string &, std::size_t, std::size_t) override;
    ~FastPositionTracker() = default;
};
