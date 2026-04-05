#pragma once

#include "types.h"

namespace eulerian_sph
{
struct WallParticle
{
    Vec2 pos = Vec2::Zero();
    Vec2 normal = Vec2::UnitX();
    Vec2 velocity_average = Vec2::Zero();
    Vec2 acceleration_average = Vec2::Zero();
    Scalar volume = 1.0;
    Scalar signed_distance = 0.0;
};

class WallParticleSet
{
  public:
    using Container = std::vector<WallParticle>;

    void reserve(const std::size_t count) { particles_.reserve(count); }
    void clear() { particles_.clear(); }
    bool empty() const { return particles_.empty(); }
    std::size_t size() const { return particles_.size(); }

    void addParticle(const WallParticle &particle) { particles_.push_back(particle); }

    WallParticle &operator[](const std::size_t index) { return particles_[index]; }
    const WallParticle &operator[](const std::size_t index) const { return particles_[index]; }

    Container &data() { return particles_; }
    const Container &data() const { return particles_; }

  private:
    Container particles_;
};
} // namespace eulerian_sph
