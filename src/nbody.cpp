#include "nbody.h"

#include <cmath>
#include <fstream>
#include <iterator>

Cartesian::Cartesian(double x, double y)
    : x(x)
    , y(y)
{
}

Cartesian & Cartesian::operator+=(const Cartesian & rhs)
{
    x += rhs.x;
    y += rhs.y;
    return *this;
}

Cartesian operator+(Cartesian const & lhs, Cartesian const & rhs)
{
    return {lhs.x + rhs.x, lhs.y + rhs.y};
}

Cartesian operator-(Cartesian const & lhs, Cartesian const & rhs)
{
    return {lhs.x - rhs.x, lhs.y - rhs.y};
}

Cartesian Cartesian::operator*(double number) const
{
    return {x * number, y * number};
}

Cartesian Cartesian::operator/(double number) const
{
    return {x / number, y / number};
}

Quadrant::Quadrant(Cartesian center, double length)
    : m_center(center)
    , m_length(length)
{
}

bool Quadrant::contains(const Cartesian & p) const
{
    return (m_center.x - m_length / 2 <= p.x) && (p.x <= m_center.x + m_length / 2) && (m_center.y - m_length / 2 <= p.y) && (p.y <= m_center.y + m_length / 2);
}

double Quadrant::length() const
{
    return m_length;
}

Quadrant Quadrant::nw() const
{
    double len = m_length / 4;
    return Quadrant({m_center.x - len, m_center.y + len}, len * 2);
}

Quadrant Quadrant::ne() const
{
    double len = m_length / 4;
    return Quadrant({m_center.x + len, m_center.y + len}, len * 2);
}

Quadrant Quadrant::sw() const
{
    double len = m_length / 4;
    return Quadrant({m_center.x - len, m_center.y - len}, len * 2);
}

Quadrant Quadrant::se() const
{
    double len = m_length / 4;
    return Quadrant({m_center.x + len, m_center.y - len}, len * 2);
}

std::ostream & operator<<(std::ostream & out, const Quadrant & b)
{
    return out << b.m_center.x << ' ' << b.m_center.y << ' ' << b.length();
}

Body::Body(const std::string & name, double weight, Cartesian coord, Cartesian velocity)
    : m_name(name)
    , m_weight(weight)
    , m_coord(coord)
    , m_velocity(velocity)
{
}

std::string Body::getName() const
{
    return m_name;
}

double Body::getWeight() const
{
    return m_weight;
}

Cartesian Body::getCoord() const
{
    return m_coord;
}

Cartesian Body::getForse() const
{
    return m_force;
}

Cartesian Body::getVelocity() const
{
    return m_velocity;
}

double Body::distance(const Body & body) const
{
    return std::hypot(m_coord.x - body.getCoord().x, m_coord.y - body.getCoord().y);
}

void Body::add_force(const Body & body)
{
    double r = distance(body);
    if (r == 0) {
        return;
    }
    double F = G * (m_weight / r) * (body.getWeight() / r);
    m_force += (body.getCoord() - m_coord) * F / r;
}

void Body::reset_force()
{
    m_force.x = 0;
    m_force.y = 0;
}

void Body::update(double delta_t)
{
    Cartesian acceleration = m_force / m_weight;
    m_velocity += acceleration * delta_t;
    m_coord += m_velocity * delta_t;
}

std::ostream & operator<<(std::ostream & out, const Body & b)
{
    return out << b.getCoord().x << ' ' << b.getCoord().y << ' ' << b.getVelocity().x << ' ' << b.getVelocity().y << ' ' << b.getWeight() << ' ' << b.getName();
}

std::istream & operator>>(std::istream & inp, Body & b)
{
    return inp >> b.m_coord.x >> b.m_coord.y >> b.m_velocity.x >> b.m_velocity.y >> b.m_weight >> b.m_name;
}

bool Body::in(const Quadrant q) const
{
    return q.contains(m_coord);
}

Body Body::plus(const Body & b) const
{
    double m = m_weight + b.getWeight();
    Cartesian coord = (m_coord * m_weight + b.getCoord() * b.getWeight()) / m;
    Cartesian velocity = (m_velocity * m_weight + b.getVelocity() * b.getWeight()) / m;
    return Body(m_name, m, coord, velocity);
}

BHTreeNode::BHTreeNode(const Quadrant & quadrant)
    : m_borders(quadrant)
{
}

void BHTreeNode::insert(const Body & b)
{
    if (!m_body && !m_nw) {
        m_body = std::make_unique<Body>(b);
        m_total_weight = b.getWeight();
        m_mass_center = b.getCoord();
    }
    else {
        if (!m_nw) {
            m_nw = std::make_unique<BHTreeNode>(m_borders.nw());
            m_ne = std::make_unique<BHTreeNode>(m_borders.ne());
            m_sw = std::make_unique<BHTreeNode>(m_borders.sw());
            m_se = std::make_unique<BHTreeNode>(m_borders.se());
        }

        m_mass_center = (m_mass_center * m_total_weight + b.getCoord() * b.getWeight()) / (m_total_weight + b.getWeight());
        m_total_weight += b.getWeight();

        if (m_body) {
            std::unique_ptr<BHTreeNode> & quad = select_quad(*m_body);
            quad->insert(*m_body);
            m_body.reset();
        }
        std::unique_ptr<BHTreeNode> & quad = select_quad(b);
        quad->insert(b);
    }
}

std::unique_ptr<BHTreeNode> & BHTreeNode::select_quad(const Body & b)
{
    if (b.in(m_nw->m_borders)) {
        return m_nw;
    }
    else if (b.in(m_ne->m_borders)) {
        return m_ne;
    }
    else if (b.in(m_sw->m_borders)) {
        return m_sw;
    }
    else {
        return m_se;
    }
}

void BHTreeNode::update_force(Body & b)
{
    if (m_body) {
        if (b.distance(*m_body) != 0) {
            b.add_force(*m_body);
        }
    }
    else {
        double dist = b.distance(Body("", m_total_weight, m_mass_center, {0, 0}));
        if (dist != 0 && m_borders.length() / dist < Theta) {
            b.add_force(Body("", m_total_weight, m_mass_center, {0, 0}));
        }
        else if (m_nw) {
            m_nw->update_force(b);
            m_ne->update_force(b);
            m_sw->update_force(b);
            m_se->update_force(b);
        }
    }
}

PositionTracker::PositionTracker(const std::string & filename)
{
    std::ifstream input(filename);
    input >> m_size;
    std::istream_iterator<Body> inp_iterator(input);
    std::istream_iterator<Body> end_input;
    std::vector<Body> bodies(inp_iterator, end_input);
    m_bodies = std::move(bodies);
    input.close();
}

BasicPositionTracker::BasicPositionTracker(const std::string & filename)
    : PositionTracker(filename)
{
}

Track BasicPositionTracker::track(const std::string & body_name, std::size_t end_time, std::size_t time_step)
{
    std::size_t ind = 0;
    while (ind < m_bodies.size() && body_name != m_bodies[ind].getName()) {
        ++ind;
    }

    Track result;
    result.push_back(m_bodies[ind].getCoord());

    for (std::size_t cur_time = 0; cur_time < end_time; cur_time += time_step) {
        for (auto & first : m_bodies) {
            for (const auto & second : m_bodies) {
                first.add_force(second);
            }
        }

        for (auto & body : m_bodies) {
            body.update(time_step);
            body.reset_force();
        }
        result.push_back(m_bodies[ind].getCoord());
    }
    return result;
}

FastPositionTracker::FastPositionTracker(const std::string & filename)
    : PositionTracker(filename)
{
}

Track FastPositionTracker::track(const std::string & body_name, std::size_t end_time, std::size_t time_step)
{
    std::size_t ind = 0;
    while (ind < m_bodies.size() && body_name != m_bodies[ind].getName()) {
        ++ind;
    }

    Track result;
    result.push_back(m_bodies[ind].getCoord());

    for (std::size_t cur_time = 0; cur_time < end_time; cur_time += time_step) {
        m_root = std::make_unique<BHTreeNode>(Quadrant({0, 0}, m_size));

        for (const auto & body : m_bodies) {
            m_root->insert(body);
        }

        for (std::size_t i = 0; i < m_bodies.size(); ++i) {
            m_root->update_force(m_bodies[i]);
            m_bodies[i].update(time_step);
            m_bodies[i].reset_force();
        }
        result.push_back(m_bodies[ind].getCoord());

        m_root.reset();
    }
    return result;
}
