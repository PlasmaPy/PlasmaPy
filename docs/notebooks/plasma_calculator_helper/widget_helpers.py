import ipywidgets as widgets
import importlib
import json
from plasmapy.particles import Particle
from astropy.constants.si import m_p, m_e
import astropy.units as units
from inspect import signature, trace


values_container = dict()
process_queue = []

class GenericWidget:
    def __init__(self, property_name,property_alias="",values_cont=values_container):
        self.property_name = property_name
        self.property_alias = property_alias or property_name
        self.widget = None
        self.values_cont = values_cont
        self.unit = None
        self.units_dropdown = None
    def set_unit(self,unit):
        self.unit = unit
    def get_widget(self):
        return self.widget
    def get_dropdown_widget(self):
        return self.units_dropdown
    def set_place_holder(self,text):
        self.widget.placeholder = text
    def create_widget(self):
        raise NotImplementedError()
    def post_creation(self):
        self.set_place_holder("")
        self.widget.observe(self.handle_change,names="value")
    def edge_case(self,value):
        pass
    def edge_case_condition(self,value):
        return False
    def try_change_value(self,value):
        self.values_cont[self.property_name]=value
    def display_error(self,value):
        if self.widget:
            self.widget.layout.border="2px solid red"
            self.widget.description = "Invalid "+self.property_alias
            self.values_cont[self.property_name] = None
    def convert_to_unit(self,change):
        if self.unit:
            return change.new*self.unit
        return change.new
    def attach_units_dropdown(self,options):
        self.units_dropdown = widgets.Dropdown(
            options=options,
            value=options[0],
            layout=widgets.Layout(width="100px")
        )
        self.units_dropdown.observe(self.handle_dropdown_change,names="value")
    def handle_dropdown_change(self,change):
        self.set_unit(self.units_dropdown.value)
        if self.property_name in self.values_cont:
            self.values_cont[self.property_name] = self.values_cont[self.property_name].value*self.unit
    def handle_change(self,change):
        value = self.convert_to_unit(change)
        if self.edge_case_condition(value):
            self.edge_case(value)
        else:
            try:
                self.try_change_value(value)
            except:
                self.display_error(value)

class FloatBox(GenericWidget):
    def __init__(self, property_name,min=-1e50,max=1e50):
        super().__init__(property_name)
        self.min = min
        self.max = max
    def create_widget(self,style={"description_width": "initial"}):
        self.widget = widgets.BoundedFloatText(name=self.property_name,min=self.min, 
            max=self.max, value=0, step=0.1,style=style)
        self.post_creation()

class CheckBox(GenericWidget):
    def __init__(self, property_name):
        super().__init__(property_name)
    def create_widget(self):
        self.widget = widgets.Checkbox(value=False)
        self.post_creation()
    def try_change_value(self,value):
        self.values_cont[self.property_name]=value

class ParticleBox(GenericWidget):
    def __init__(self, property_name,property_alias=None):
        super().__init__(property_name,property_alias=property_alias)

    def edge_case_condition(self, value):
        return value is None or value == ""
    def edge_case(self, value):
        self.values_cont[self.property_name] = None
        self.widget.description = ""
        self.widget.layout.border = "1px solid black"
    def try_change_value(self, value):
        particle = Particle(value)
        self.values_cont[self.property_name] = particle
        self.widget.layout.border="1px solid black"
        self.widget.description = ""

    def create_widget(self,style={"description_width": "initial"}):
        self.widget = widgets.Text(style=style)
        self.post_creation()

class IonBox(ParticleBox):
    def __init__(self, property_name,property_alias=None):
        super().__init__(property_name,property_alias=property_alias)
    def try_change_value(self, value):
        ion = Particle(value)
        if ion.is_ion:
            self.values_cont[self.property_name] = ion
            self.widget.layout.border="1px solid black"
            self.widget.description = ""
        else:
            raise ValueError(f"{ion} is not an ion")

class FunctionInfo:
    def __init__(self, module_name, function_name, values_cont=values_container):
        self.module = module_name
        self.fname = function_name
        self.fattr = getattr(importlib.import_module(module_name),function_name)
        self.values_cont = values_cont
        self.spec_combo = None
        self.sig = list(signature(self.fattr).parameters.keys())
        self.output_widget = widgets.Output()
        self.output_widget.layout.margin = "10px 10px 10px 10px"
        self.output_widget.layout.padding = "10px 10px 10px 10px"
    def add_combo(self,spec_combo):
        if not self.spec_combo:
            self.spec_combo = []
        self.spec_combo.append(spec_combo)
    def get_output_widget(self):
        return self.output_widget
    def produce_arg(self, spec):
        args_dict = dict()
        for arg in spec:
            if arg in self.values_cont and self.values_cont[arg] is not None:
                args_dict[arg] = self.values_cont[arg]

        return args_dict

    def error_message(self,spec):
        print(colored(0,0,0,"["),end="")
        for arg in spec:
            if arg in self.values_cont and self.values_cont[arg] is not None:
                print(colored(0,128,0,arg+":present,"),end="")
            else:
                print(colored(255,0,0,arg+":missing,"),end="")
        print(colored(0,0,0,"]"))
    def process(self):
        self.output_widget.clear_output()
        args_dict = dict()
        if self.spec_combo:
           for spec in self.spec_combo:
               args_dict = self.produce_arg(spec)
               
               if len(args_dict) == len(spec):
                   break
        else:
            args_dict = self.produce_arg(self.sig)
        with self.output_widget:
            try:
                self.output_widget.layout.border="0px"
                print(" : "+str(self.fattr(**args_dict)))

            except Exception as e:
                self.output_widget.layout.border="1px solid red"
                print(e)
                print(" : could not be computed one or more parameter is missing - check below for missing parameters")
                if self.spec_combo:
                    for spec in self.spec_combo:
                        self.error_message(spec)
                else:
                    self.error_message(self.sig)



def create_label(label,color="black"):
    return widgets.HTML(f"<h3 style='margin:0px;color:{color}'>{label}<h3>")

def create_button():
    button = widgets.Button(description="Calculate Properties",button_style="info")
    return button

def handle_button_click(event):
    print("clicked")
    for fn in process_queue:
        fn.process()

def handle_clear_click(event):
    for key in values_container:
        values_container[key] = None
    for fn in process_queue:
        fn.output_widget.clear_output()
calculate_button = widgets.Button(description="Calculate Properties",button_style="info")
calculate_button.on_click(handle_button_click)

clear_button = widgets.Button(description="Clear Values",button_style="danger")
clear_button.on_click(handle_clear_click)

def create_widget(widget_type,**kwargs):
    unit = None
    placeholder = None
    opts=None
    if 'unit' in kwargs:
        unit = kwargs['unit']
        del kwargs['unit']
    if 'placeholder' in kwargs:
        placeholder = kwargs['placeholder']
        del kwargs['placeholder']
    if 'opts' in kwargs:
        opts = kwargs['opts']
        del kwargs['opts']
    widget_element = widget_type(**kwargs)
    widget_element.create_widget()
    if unit:
        widget_element.set_unit(unit)
    if placeholder:
        widget_element.set_place_holder(placeholder)
    if opts:
        widget_element.attach_units_dropdown(opts)
        return [widget_element.get_widget(),widget_element.get_dropdown_widget()]
    return widget_element.get_widget()




def colored(r, g, b, text):
    return "\033[38;2;{};{};{}m{} \033[38;2;255;255;255m".format(r, g, b, text)
