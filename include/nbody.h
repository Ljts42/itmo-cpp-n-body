#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string_view>
#include <unordered_map>
#include <vector>

const double G = 6.67e-11;
const double Theta = 0.5;

struct Cartesian
{
    double x = 0;
    double y = 0;

    Cartesian() = default;
    Cartesian(const Cartesian &) = default;
    Cartesian(double x, double y)
        : x(x)
        , y(y)
    {
    }

    Cartesian & operator+=(const Cartesian & rhs)
    {
        this->x += rhs.x;
        this->y += rhs.y;
        return *this;
    }
};

// Quadrant representation, required for Problem 2
class Quadrant
{
    Cartesian m_center;
    double m_length;

public:
    // Create quadrant with center (x, y) and size 'lenth'
    Quadrant(Cartesian, double);
    Quadrant(double, double, double);

    // Test if point (x, y) is in the quadrant
    bool contains(const Cartesian &) const;
    double length() const;

    // The four methods below construct new Quadrant representing sub-quadrant of the invoking quadrant
    Quadrant nw();
    Quadrant ne();
    Quadrant sw();
    Quadrant se();

    friend std::ostream & operator<<(std::ostream &, const Quadrant &);
};

// Single body representation, required for Problem 1 and Problem 2
class Body
{
    double m_weight = -1;
    Cartesian m_coord;
    Cartesian m_velocity;
    Cartesian m_force{0, 0};

public:
    double getWeight() const { return m_weight; }
    double getX() const { return m_coord.x; }
    double getY() const { return m_coord.y; }
    Cartesian getCoord() const { return m_coord; }
    Cartesian getForse() const { return m_force; }
    Cartesian getSpeed() const { return m_velocity; }

    Body() = default;
    Body(double, Cartesian, Cartesian);

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
    Body plus(const Body &);
};

// Burnes-Hut tree representation, required for Problem 2
class BHTreeNode
{
    double weight;
    Cartesian mass_center;
    Quadrant borders;
    std::shared_ptr<Body> body = nullptr;

    std::shared_ptr<BHTreeNode> northWest = nullptr;
    std::shared_ptr<BHTreeNode> northEast = nullptr;
    std::shared_ptr<BHTreeNode> southWest = nullptr;
    std::shared_ptr<BHTreeNode> southEast = nullptr;

public:
    BHTreeNode() = default;
    BHTreeNode(Quadrant qua)
        : borders(qua)
    {
    }

    void insert(std::shared_ptr<Body>);
    // Update net acting force-on 'b'
    void update_force(Body &);

    ~BHTreeNode()
    {
        clear(northWest);
        clear(northEast);
        clear(southWest);
        clear(southEast);
    }

    void clear(std::shared_ptr<BHTreeNode> node)
    {
        if (node) {
            clear(node->northWest);
            clear(node->northEast);
            clear(node->southWest);
            clear(node->southEast);
            node.reset();
        }
    }

private:
    bool hasBody() const
    {
        return body != nullptr;
    }
};

using Track = std::vector<Cartesian>;

class PositionTracker
{
protected:
    std::unordered_map<std::string, std::shared_ptr<Body>> bodies;
    double size;
    PositionTracker(const std::string &);

public:
    virtual Track track(const std::string &, std::size_t, std::size_t) = 0;
};

class BasicPositionTracker : public PositionTracker
{
public:
    BasicPositionTracker(const std::string &);
    Track track(const std::string &, std::size_t, std::size_t) override;
};

class FastPositionTracker : public PositionTracker
{
    BHTreeNode * root;

public:
    FastPositionTracker(const std::string &);
    Track track(const std::string &, std::size_t, std::size_t) override;
};
