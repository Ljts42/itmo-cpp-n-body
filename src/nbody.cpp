#include "nbody.h"

Body::Body(double weight, Cartesian coord, Cartesian velocity)
    : m_weight(weight)
    , m_coord(coord)
    , m_velocity(velocity)
{
}

double Body::distance(const Body & b) const
{
    double dx = std::abs(m_coord.x - b.m_coord.x);
    double dy = std::abs(m_coord.y - b.m_coord.y);
    return std::hypot(dx, dy);
}

void Body::add_force(const Body & b)
{
    double r = distance(b);
    if (r == 0) {
        return;
    }
    double F = G * (m_weight / r) * (b.m_weight / r);
    m_force += Cartesian(F * ((b.m_coord.x - m_coord.x) / r),
                         F * ((b.m_coord.y - m_coord.y) / r));
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
    return out << b.m_coord.x << ' ' << b.m_coord.y << ' ' << b.m_velocity.x << ' ' << b.m_velocity.y << ' ' << b.m_weight;
}

std::istream & operator>>(std::istream & inp, Body & b)
{
    return inp >> b.m_coord.x >> b.m_coord.y >> b.m_velocity.x >> b.m_velocity.y >> b.m_weight;
}

bool Body::in(const Quadrant q) const
{
    return q.contains(m_coord);
}

Body Body::plus(const Body & b) const
{
    double m = m_weight + b.m_weight;

    double x = m_coord.x * (m_weight / m) + b.m_coord.x * (b.m_weight / m);
    double y = m_coord.y * (m_weight / m) + b.m_coord.y * (b.m_weight / m);

    double v_x = m_velocity.x * (m_weight / m) + b.m_velocity.x * (b.m_weight / m);
    double v_y = m_velocity.y * (m_weight / m) + b.m_velocity.y * (b.m_weight / m);
    return Body(m, Cartesian(x, y), Cartesian(v_x, v_y));
}

Quadrant::Quadrant(Cartesian center, double length)
    : m_center(center)
    , m_length(length)
{
}

Quadrant::Quadrant(double x, double y, double length)
    : m_center(x, y)
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
    return Quadrant(m_center.x - len, m_center.y + len, len);
}

Quadrant Quadrant::ne() const
{
    double len = m_length / 2;
    return Quadrant(m_center.x + len, m_center.y + len, len);
}

Quadrant Quadrant::sw() const
{
    double len = m_length / 2;
    return Quadrant(m_center.x - len, m_center.y - len, len);
}

Quadrant Quadrant::se() const
{
    double len = m_length / 2;
    return Quadrant(m_center.x + len, m_center.y - len, len);
}

std::ostream & operator<<(std::ostream & out, const Quadrant & b)
{
    return out << b.m_center.x << ' ' << b.m_center.y << ' ' << b.m_length;
}

void BHTreeNode::insert(std::shared_ptr<Body> & b)
{
    if (!hasBody() && northWest == nullptr) {
        body = b;
        weight = b->getWeight();
        mass_center = b->getCoord();
    }
    else {
        if (northWest == nullptr) {
            northWest = std::make_shared<BHTreeNode>(borders.nw());
            northEast = std::make_shared<BHTreeNode>(borders.ne());
            southWest = std::make_shared<BHTreeNode>(borders.sw());
            southEast = std::make_shared<BHTreeNode>(borders.se());
        }

        double x = mass_center.x * weight + b->getX() * b->getWeight();
        double y = mass_center.y * weight + b->getY() * b->getWeight();
        weight += b->getWeight();
        mass_center = Cartesian(x / weight, y / weight);

        if (hasBody()) {
            if (body->in(northWest->borders) && !b->in(northWest->borders)) {
                northWest->insert(body);
            }
            else if (body->in(northEast->borders) && !b->in(northWest->borders)) {
                northEast->insert(body);
            }
            else if (body->in(southWest->borders) && !b->in(northWest->borders)) {
                southWest->insert(body);
            }
            else if (body->in(southEast->borders) && !b->in(southEast->borders)) {
                southEast->insert(body);
            }
            else if (body->in(southEast->borders)) {
                southEast->insert(body);
            }
            else if (body->in(southWest->borders)) {
                southWest->insert(body);
            }
            else if (body->in(northEast->borders)) {
                northEast->insert(body);
            }
            else {
                northWest->insert(body);
            }
            body = nullptr;
        }
        if (b->in(northWest->borders)) {
            northWest->insert(b);
        }
        else if (b->in(northEast->borders)) {
            northEast->insert(b);
        }
        else if (b->in(southWest->borders)) {
            southWest->insert(b);
        }
        else {
            southEast->insert(b);
        }
    }
}

void BHTreeNode::update_force(Body & b)
{
    if (hasBody()) {
        if (b.distance(*body) != 0) {
            b.add_force(*body);
        }
    }
    else {
        double coef = borders.length() / b.distance(Body(weight, mass_center, {0, 0}));
        if (coef < Theta) {
            b.add_force(Body(weight, mass_center, {0, 0}));
        }
        else if (northWest != nullptr) {
            northWest->update_force(b);
            northEast->update_force(b);
            southWest->update_force(b);
            southEast->update_force(b);
        }
    }
}

PositionTracker::PositionTracker(const std::string & filename)
{
    std::ifstream input(filename);
    double size;
    std::string c;
    input >> size;
    std::getline(input, c);
    Body b;
    std::string name;
    while (input >> b >> name) {
        std::getline(input, c);
        bodies[name] = std::make_shared<Body>(b);
    }
}

BasicPositionTracker::BasicPositionTracker(const std::string & filename)
    : PositionTracker(filename)
{
}

Track BasicPositionTracker::track(const std::string & body_name, std::size_t end_time, std::size_t time_step)
{
    Track result;
    result.push_back(bodies[body_name]->getCoord());
    for (std::size_t cur_time = 0; cur_time < end_time; cur_time += time_step) {
        for (auto & first : bodies) {
            for (auto & second : bodies) {
                first.second->add_force(*second.second);
            }
        }
        for (auto & body : bodies) {
            body.second->update(time_step);
            body.second->reset_force();
        }
        result.push_back(bodies[body_name]->getCoord());
    }
    return result;
}

FastPositionTracker::FastPositionTracker(const std::string & filename)
    : PositionTracker(filename)
{
}

Track FastPositionTracker::track(const std::string & body_name, std::size_t end_time, std::size_t time_step)
{
    Track result(0);
    result.push_back(bodies[body_name]->getCoord());
    for (std::size_t cur_time = 0; cur_time < end_time + 54; cur_time += time_step) {
        root = new BHTreeNode(Quadrant(0, 0, size));
        for (auto & body : bodies) {
            root->insert(body.second);
        }
        for (auto & body : bodies) {
            root->update_force(*body.second);
        }
        for (auto & body : bodies) {
            body.second->update(time_step);
            body.second->reset_force();
        }
        result.push_back(bodies[body_name]->getCoord());

        delete root;
    }
    return result;
}
