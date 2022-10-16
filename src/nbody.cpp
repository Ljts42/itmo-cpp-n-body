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
    double len = m_length / 2;
    return Quadrant({m_center.x - len, m_center.y + len}, len);
}

Quadrant Quadrant::ne() const
{
    double len = m_length / 2;
    return Quadrant({m_center.x + len, m_center.y + len}, len);
}

Quadrant Quadrant::sw() const
{
    double len = m_length / 2;
    return Quadrant({m_center.x - len, m_center.y - len}, len);
}

Quadrant Quadrant::se() const
{
    double len = m_length / 2;
    return Quadrant({m_center.x + len, m_center.y - len}, len);
}

std::ostream & operator<<(std::ostream & out, const Quadrant & b)
{
    return out << b.m_center.x << ' ' << b.m_center.y << ' ' << b.m_length;
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

Cartesian Body::getSpeed() const
{
    return m_velocity;
}

double Body::distance(const Body & body) const
{
    return std::hypot(m_coord.x - body.m_coord.x, m_coord.y - body.m_coord.y);
}

void Body::add_force(const Body & body)
{
    double r = distance(body);
    if (r == 0) {
        return;
    }
    double F = G * (m_weight / r) * (body.m_weight / r);
    m_force += Cartesian(F * (body.m_coord.x - m_coord.x) / r,
                         F * (body.m_coord.y - m_coord.y) / r);
}

void Body::reset_force()
{
    m_force.x = 0;
    m_force.y = 0;
}

void Body::update(double delta_t)
{
    double a_x = m_force.x / m_weight;
    double a_y = m_force.y / m_weight;
    m_velocity += Cartesian(a_x * delta_t, a_y * delta_t);
    m_coord += Cartesian(m_velocity.x * delta_t, m_velocity.y * delta_t);
}

std::ostream & operator<<(std::ostream & out, const Body & b)
{
    return out << b.m_coord.x << ' ' << b.m_coord.y << ' ' << b.m_velocity.x << ' ' << b.m_velocity.y << ' ' << b.m_weight << ' ' << b.m_name;
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
    double m = m_weight + b.m_weight;
    double x = (m_coord.x * m_weight + b.m_coord.x * b.m_weight) / m;
    double y = (m_coord.y * m_weight + b.m_coord.y * b.m_weight) / m;

    double v_x = m_velocity.x * (m_weight / m) + b.m_velocity.x * (b.m_weight / m);
    double v_y = m_velocity.y * (m_weight / m) + b.m_velocity.y * (b.m_weight / m);
    return Body(m_name, m, Cartesian(x, y), Cartesian(v_x, v_y));
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

        double x = m_mass_center.x * m_total_weight + b.getCoord().x * b.getWeight();
        double y = m_mass_center.y * m_total_weight + b.getCoord().y * b.getWeight();
        m_total_weight += b.getWeight();
        m_mass_center = Cartesian(x / m_total_weight, y / m_total_weight);

        if (m_body) {
            if (m_body->in(m_nw->m_borders) && !b.in(m_nw->m_borders)) {
                m_nw->insert(*m_body);
            }
            else if (m_body->in(m_ne->m_borders) && !b.in(m_nw->m_borders)) {
                m_ne->insert(*m_body);
            }
            else if (m_body->in(m_sw->m_borders) && !b.in(m_nw->m_borders)) {
                m_sw->insert(*m_body);
            }
            else if (m_body->in(m_se->m_borders) && !b.in(m_se->m_borders)) {
                m_se->insert(*m_body);
            }
            else if (m_body->in(m_se->m_borders)) {
                m_se->insert(*m_body);
            }
            else if (m_body->in(m_sw->m_borders)) {
                m_sw->insert(*m_body);
            }
            else if (m_body->in(m_ne->m_borders)) {
                m_ne->insert(*m_body);
            }
            else {
                m_nw->insert(*m_body);
            }
            m_body = nullptr;
        }
        if (b.in(m_nw->m_borders)) {
            m_nw->insert(b);
        }
        else if (b.in(m_ne->m_borders)) {
            m_ne->insert(b);
        }
        else if (b.in(m_sw->m_borders)) {
            m_sw->insert(b);
        }
        else {
            m_se->insert(b);
        }
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
        m_root = new BHTreeNode(Quadrant({0, 0}, m_size));

        for (const auto & body : m_bodies) {
            m_root->insert(body);
        }
        for (std::size_t i = 0; i < m_bodies.size(); ++i) {
            m_root->update_force(m_bodies[i]);
            m_bodies[i].update(time_step);
            m_bodies[i].reset_force();
        }
        result.push_back(m_bodies[ind].getCoord());
        delete m_root;
        m_root = nullptr;
    }
    return result;
}
