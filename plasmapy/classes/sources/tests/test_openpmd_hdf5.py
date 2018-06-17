import pytest
import plasmapy
from astropy import units as u
import os

@pytest.fixture(scope="module")
def openPMD2DPlasma():
    """This is me not knowing how to handle the new Plasma class with data!

    Downloaded from https://github.com/openPMD/openPMD-example-datasets/blob/draft/example-2d.tar.gz
    per the Creative Commons Zero v1.0 Universal license"""
    return plasmapy.classes.Plasma(
        hdf5=os.path.join(
            os.path.dirname(__file__) or '.',
            "data",
            "data00000255.h5")
        )

def test_has_electric_field_with_units(openPMD2DPlasma):
    "It would be nice to let the fields have the correct units from the get go"
    assert openPMD2DPlasma.electric_field.to(u.V / u.m)

def test_has_charge_density_with_units(openPMD2DPlasma):
    assert openPMD2DPlasma.charge_density.to(u.C/u.m**3)  # unless it's some
    # 2D charge density

def test_correct_shape_electric_field(openPMD2DPlasma):
    assert openPMD2DPlasma.electric_field.shape == (3, 51, 201)
    # IIRC this is how we defined it in the old Plasma class but I can be
    # convinced otherwise if another way makes more sense

def test_x(openPMD2DPlasma):
    assert openPMD2DPlasma.x.shape == (51, 201)
    assert openPMD2DPlasma.x.to(u.m)

def test_y(openPMD2DPlasma):
    """"I'm not sure this should happen but the 2D dataset has an x-z grid only so
    we will have to discuss this"""
    with pytest.raises(AttributeError):
        print(openPMD2DPlasma.y)

def test_z(openPMD2DPlasma):
    assert openPMD2DPlasma.z.shape == (51, 201)
    assert openPMD2DPlasma.z.to(u.m)


