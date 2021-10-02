import ipywidgets as widgets
import importlib
from plasmapy.particles import Particle
from astropy.constants.si import m_p, m_e
import astropy.units as units
from inspect import signature, trace

BLACK_COLOR = (0,0,0)
DARK_RED_COLOR = (255,0,0)
LIGHT_GREEN_COLOR = (0,128,0)

ERROR_STYLE = "2px solid red"
EQUAL_SPACING_CONFIG = "10px 10px 10px 10px"

values_container = dict()
_process_queue = []

class _GenericWidget:
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
            self.widget.layout.border = ERROR_STYLE
            self.widget.description = f"Invalid ${self.property_alias}"
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

class _FloatBox(_GenericWidget):
    def __init__(self, property_name,min=-1e50,max=1e50):
        super().__init__(property_name)
        self.min = min
        self.max = max
    def create_widget(self,style = {"description_width": "initial"}):
        self.widget = widgets.BoundedFloatText(name=self.property_name,min=self.min, 
            max=self.max, value=0, step=0.1,style=style)
        self.post_creation()

class _CheckBox(_GenericWidget):
    def __init__(self, property_name):
        super().__init__(property_name)
    def create_widget(self):
        self.widget = widgets.Checkbox(value=False)
        self.post_creation()
    def try_change_value(self,value):
        self.values_cont[self.property_name]=value

class _ParticleBox(_GenericWidget):
    def __init__(self, property_name,property_alias=None):
        super().__init__(property_name,property_alias=property_alias)

    def edge_case_condition(self, value):
        return value is None or value == ""
    def edge_case(self, value):
        self.values_cont[self.property_name] = None
        self.widget.description = ""
        self.widget.layout.border = ""
    def try_change_value(self, value):
        particle = Particle(value)
        self.values_cont[self.property_name] = particle
        self.widget.layout.border=""
        self.widget.description = ""

    def create_widget(self,style = {"description_width": "initial"}):
        self.widget = widgets.Text(style=style)
        self.post_creation()

class _IonBox(_ParticleBox):
    def __init__(self, property_name,property_alias=None):
        super().__init__(property_name,property_alias=property_alias)
    def try_change_value(self, value):
        ion = Particle(value)
        if ion.is_ion:
            self.values_cont[self.property_name] = ion
            self.widget.layout.border=""
            self.widget.description = ""
        else:
            raise ValueError(f"{ion} is not an ion")

class _FunctionInfo:
    def __init__(self, module_name, function_name, values_cont=values_container):
        self.module = module_name
        self.fname = function_name
        self.fattr = getattr(importlib.import_module(module_name),function_name)
        self.values_cont = values_cont
        self.spec_combo = None
        self.sig = list(signature(self.fattr).parameters.keys())
        self.output_widget = widgets.Output()
        self.output_widget.layout.margin = EQUAL_SPACING_CONFIG
        self.output_widget.layout.padding = EQUAL_SPACING_CONFIG
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
        print(_colored_text(BLACK_COLOR, "["),end="")
        for arg in spec:
            if arg in self.values_cont and self.values_cont[arg] is not None:
                print(_colored_text(LIGHT_GREEN_COLOR, arg+":present,"),end="")
            else:
                print(_colored_text(DARK_RED_COLOR, arg+":missing,"),end="")
        print(_colored_text(BLACK_COLOR, "]"))
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
                self.output_widget.layout.border = "0px"
                print(" : "+str(self.fattr(**args_dict)))

            except Exception as e:
                self.output_widget.layout.border = ERROR_STYLE
                print(e)
                print(" : could not be computed one or more parameter is missing - check below for missing parameters")
                if self.spec_combo:
                    for spec in self.spec_combo:
                        self.error_message(spec)
                else:
                    self.error_message(self.sig)


def _create_label(label,color="black"):
    return widgets.HTML(f"<h3 style='margin:0px;color:{color}'>{label}<h3>")

def _create_button():
    button = widgets.Button(description="Calculate Properties",button_style="info")
    return button

def _handle_button_click(event):
    for fn in _process_queue:
        fn.process()

def _handle_clear_click(event):
    for fn in _process_queue:
        fn.output_widget.clear_output()

def _colored_text(color, text):
    return "\033[38;2;{};{};{}m{} \033[38;2;255;255;255m".format(color[0], color[1], color[2], text)

_calculate_button = widgets.Button(description="Calculate Properties",button_style="info")
_calculate_button.on_click(_handle_button_click)

_clear_button = widgets.Button(description="Clear Output",button_style="danger")
_clear_button.on_click(_handle_clear_click)

def _create_widget(widget_type,**kwargs):
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
        widgets = [widget_element.get_widget(),widget_element.get_dropdown_widget()]   
    else:
        widgets = widget_element.get_widget()
    return widgets
