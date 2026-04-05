#pragma once

#include "types.h"

namespace eulerian_sph
{
struct ParticleState
{
    Vec2 pos = Vec2::Zero();
    Vec2 vel = Vec2::Zero();
    Vec2 mom = Vec2::Zero();
    Vec2 normal = Vec2::UnitX();
    Scalar rho = 1.0;
    Scalar p = 0.0;
    Scalar mass = 0.0;
    Scalar volume = 1.0;
    Scalar energy = 0.0;
    Scalar signed_distance = 0.0;
    Scalar position_divergence = 0.0;
    int indicator = 0;
    int smeared_surface = 0;
    int boundary_type = 0;
};

class ParticleSet
{
  public:
    using Container = std::vector<ParticleState>;
    using iterator = Container::iterator;
    using const_iterator = Container::const_iterator;

    void reserve(const std::size_t count) { particles_.reserve(count); }
    void clear() { particles_.clear(); }
    bool empty() const { return particles_.empty(); }
    std::size_t size() const { return particles_.size(); }

    void addParticle(const ParticleState &particle) { particles_.push_back(particle); }

    ParticleState &operator[](const std::size_t index) { return particles_[index]; }
    const ParticleState &operator[](const std::size_t index) const { return particles_[index]; }

    iterator begin() { return particles_.begin(); }
    iterator end() { return particles_.end(); }
    const_iterator begin() const { return particles_.begin(); }
    const_iterator end() const { return particles_.end(); }

    Container &data() { return particles_; }
    const Container &data() const { return particles_; }

  private:
    Container particles_;
};

inline void synchronizeConservativeVariables(ParticleState &particle)
{
    particle.mass = particle.rho * particle.volume;
    particle.mom = particle.mass * particle.vel;
}
} // namespace eulerian_sph
