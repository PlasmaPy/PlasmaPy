from IPython.display import display
from .widget_helpers import *

grid_data = [
    [
        create_label("Parameter",color="#00BFD8"),
        create_label("Value",color="#00BFD8"),
        create_label("Unit",color="#00BFD8")
    ],
    [
        create_label("B - Magnetic Field Magnitude:"),
        *create_widget(FloatBox,property_name="B",unit=units.T,opts=[units.T,units.G,units.uG])
    ],
    [
        create_label("Particle:"),
        create_widget(ParticleBox,property_name="particle",property_alias="particle_type",placeholder="Enter particle e.g. neutron,He")
    ],
    [
        create_label("Ion:"),
        create_widget(IonBox,property_name="ion",property_alias="ion_type",placeholder="Enter ion e.g. He 2+")
    ],
    [
        create_label("Convert to Hertz:"),
        create_widget(CheckBox,property_name="to_hz")
    ],
    [
        create_label("_"*20,color="#00BFD8"),
        create_label("Density Number"),
        create_label("_"*20,color="#00BFD8")
    ],
    [
        create_label("n - Standard Density Number:"),
        *create_widget(FloatBox,property_name="n",unit=units.m**-3,
            opts=[units.m**-3,units.cm**-3,units.mm**-3])
    ],
    [
        create_label("n<sub>e</sub> - Electron Density Number:"),
        *create_widget(FloatBox,property_name="n_e",unit=units.m**-3,
            opts=[units.m**-3,units.cm**-3,units.mm**-3])
    ],
    [
        create_label("n<sub>i</sub> - Ion Density Number:"),
        *create_widget(FloatBox,property_name="n_i",unit=units.m**-3,
            opts=[units.m**-3,units.cm**-3,units.mm**-3])
    ],
    [
        create_label("_"*20,color="#00BFD8"),
        create_label("Temperature"),
        create_label("_"*20,color="#00BFD8")
    ],
    [
        create_label("T - Standard Temperature:"),
        create_widget(FloatBox,property_name="T",min=0,unit=units.K),
        create_label("K",color="#a9a9a9")
    ],
    [
        create_label("T<sub>e</sub> - Electron Temperature:"),
        create_widget(FloatBox,property_name="T_e",min=0,unit=units.K),
        create_label("K",color="#a9a9a9")
    ],
    [
        create_label("T<sub>i</sub> - Ion Temperature:"),
        create_widget(FloatBox,property_name="T_i",min=0,unit=units.K),
        create_label("K",color="#a9a9a9")
    ],
    
]




def create_interactive_layout():
    ## grid config
    grid = widgets.GridspecLayout(15, 3)
    grid.layout.margin="10px"

    for i,row in enumerate(grid_data):
        for j,cell in enumerate(row):
            grid[i,j] = cell
    
    grid[-1,0] = calculate_button
    return grid

def create_output_layout():
    app = widgets.Tab()
    children = []

    with open('plasma_calculator_helper/properties_metadata.json') as f:
        data = json.load(f)

    for i,title in enumerate(data):
        grid_layout = widgets.GridspecLayout(10,2,width="100%")
        for j,prop in enumerate(data[title]):
            fn = FunctionInfo(prop["module_name"],prop["function_name"])
            if "spec_combo" in prop:
                for spec_combo in prop["spec_combo"]:
                    fn.add_combo(spec_combo)
            grid_layout[j,0] = create_label(prop["function_name"]+":")
            grid_layout[j,1] = fn.get_output_widget()
            process_queue.append(fn)
        children.append(grid_layout)
        app.set_title(i,title)
    app.children = children

    return app

