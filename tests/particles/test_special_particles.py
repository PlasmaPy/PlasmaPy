from plasmapy.particles._special_particles import data_about_special_particles


def test_particle_antiparticle_pairs(particle_antiparticle_pair) -> None:
    """Test that particles and antiparticles have the same or exact
    opposite properties in the _Particles dictionary."""
    particle, antiparticle = (p.symbol for p in particle_antiparticle_pair)

    assert not data_about_special_particles[particle][
        "antimatter"
    ], f"{particle} is incorrectly marked as antimatter."

    assert data_about_special_particles[antiparticle][
        "antimatter"
    ], f"{antiparticle} is incorrectly marked as matter."

    identical_keys = ["half-life", "spin"]

    if "nu" not in particle:
        identical_keys.append("mass")

    if particle in ("e-", "mu-", "tau-") or "nu" in particle:
        identical_keys.append("generation")

    opposite_keys = ["charge number", "lepton number", "baryon number"]

    for key in identical_keys:
        assert (
            data_about_special_particles[particle][key]
            == data_about_special_particles[antiparticle][key]
        ), f"{particle} and {antiparticle} do not have identical {key}."

    for key in opposite_keys:
        assert (
            data_about_special_particles[particle][key]
            == -data_about_special_particles[antiparticle][key]
        ), f"{particle} and {antiparticle} do not have exact opposite {key}."

    if particle not in ("e-", "n"):
        assert data_about_special_particles[particle][
            "name"
        ] == data_about_special_particles[antiparticle]["name"].replace(
            "anti", ""
        ), f"{particle} and {antiparticle} do not have same name except for 'anti'."


required_keys = [
    "name",
    "spin",
    "class",
    "lepton number",
    "baryon number",
    "charge number",
    "half-life",
    "mass",
    "antimatter",
]


def test_Particles_required_keys(particle):
    r"""Test that required keys are present for all particles."""

    missing_keys = []

    for key in required_keys:
        try:
            data_about_special_particles[particle.symbol][key]
        except KeyError:
            missing_keys.append(key)

    if missing_keys:
        raise KeyError(
            f"_Particles[{particle!r}] is missing the following "
            f"keys: {missing_keys}"
        )
