"""
Collection of private functions to load properties and construct widgets
"""
__all__ = []

import astropy.units as units
import json

from ipywidgets import widgets
from pathlib import Path

from plasmapy.utils.calculator.widget_helpers import (
    _calculate_button,
    _CheckBox,
    _clear_button,
    _create_label,
    _create_widget,
    _FloatBox,
    _FunctionInfo,
    _IonBox,
    _ParticleBox,
    _process_queue,
)

LIGHT_BLUE = "#00BFD8"
"""
LIGHT_BLUE: Constant for light blue color 00BFD8
"""
LIGHT_GRAY = "#A9A9A9"
"""
LIGHT_GRAY: Constant for light gray color A9A9A9
"""

grid_data = [
    [
        _create_label("Parameter", color=LIGHT_BLUE),
        _create_label("Value", color=LIGHT_BLUE),
        _create_label("Unit", color=LIGHT_BLUE),
    ],
    [
        _create_label("B - Magnetic Field Magnitude:"),
        *_create_widget(
            _FloatBox,
            property_name="B",
            unit=units.T,
            opts=[units.T, units.G, units.uG],
        ),
    ],
    [
        _create_label("Particle:"),
        _create_widget(
            _ParticleBox,
            property_name="particle",
            property_alias="particle_type",
            placeholder="Enter particle e.g. neutron,He",
        ),
    ],
    [
        _create_label("Ion:"),
        _create_widget(
            _IonBox,
            property_name="ion",
            property_alias="ion_type",
            placeholder="Enter ion e.g. He 2+",
        ),
    ],
    [
        _create_label("Convert to Hertz:"),
        _create_widget(_CheckBox, property_name="to_hz"),
    ],
    [
        _create_label("_" * 20, color=LIGHT_BLUE),
        _create_label("Density Number"),
        _create_label("_" * 20, color=LIGHT_BLUE),
    ],
    [
        _create_label("n - Standard Density Number:"),
        *_create_widget(
            _FloatBox,
            property_name="n",
            unit=units.m**-3,
            opts=[units.m**-3, units.cm**-3, units.mm**-3],
        ),
    ],
    [
        _create_label("n<sub>e</sub> - Electron Density Number:"),
        *_create_widget(
            _FloatBox,
            property_name="n_e",
            unit=units.m**-3,
            opts=[units.m**-3, units.cm**-3, units.mm**-3],
        ),
    ],
    [
        _create_label("n<sub>i</sub> - Ion Number Density:"),
        *_create_widget(
            _FloatBox,
            property_name="n_i",
            unit=units.m**-3,
            opts=[units.m**-3, units.cm**-3, units.mm**-3],
        ),
    ],
    [
        _create_label("_" * 20, color=LIGHT_BLUE),
        _create_label("Temperature"),
        _create_label("_" * 20, color=LIGHT_BLUE),
    ],
    [
        _create_label("T - Standard Temperature:"),
        _create_widget(_FloatBox, property_name="T", min=0, unit=units.K),
        _create_label("K", color=LIGHT_GRAY),
    ],
    [
        _create_label("T<sub>e</sub> - Electron Temperature:"),
        _create_widget(_FloatBox, property_name="T_e", min=0, unit=units.K),
        _create_label("K", color=LIGHT_GRAY),
    ],
    [
        _create_label("T<sub>i</sub> - Ion Temperature:"),
        _create_widget(_FloatBox, property_name="T_i", min=0, unit=units.K),
        _create_label("K", color=LIGHT_GRAY),
    ],
]
"""
grid_data: Contains widgets layout for input parameters
"""


def _create_interactive_layout():
    """
    Interactive grid layout for input parameters populated in grid_data
    """
    grid = widgets.GridspecLayout(18, 3)
    grid.layout.margin = "10px"

    for i, row in enumerate(grid_data):
        for j, cell in enumerate(row):
            grid[i, j] = cell

    grid[-1, 0] = _calculate_button
    grid[-1, 1] = _clear_button
    return grid


def _create_output_layout():
    """
    Tab layout for output parameters, populated from properties_metadata.json
    """
    app = widgets.Tab()
    children = []

    properties_metadata = Path("properties_metadata.json")

    with properties_metadata.open() as f:
        data = json.load(f)

    for i, title in enumerate(data):
        grid_layout = widgets.GridspecLayout(10, 2, width="100%")
        for j, prop in enumerate(data[title]):
            fn = _FunctionInfo(prop["module_name"], prop["function_name"])
            if "spec_combo" in prop:
                for spec_combo in prop["spec_combo"]:
                    fn.add_combo(spec_combo)
            grid_layout[j, 0] = _create_label(prop["function_name"] + ":")
            grid_layout[j, 1] = fn.get_output_widget()
            _process_queue.append(fn)
        children.append(grid_layout)
        app.set_title(i, title)
    app.children = children

    return app
