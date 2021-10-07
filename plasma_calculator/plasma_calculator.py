
## Script file to launch the notebook as a standalone clean app

#voila docs/notebooks/plasma_calculator.ipynb
import sys
import os
import argparse
import plasmapy

def main():
    parser = argparse.ArgumentParser(description='Plasma calculator')
    parser.add_argument('--port', type=int, default=8866, help='Port to run the notebook')
    parser.add_argument('--no-gui', action='store_true', help='Do not open the GUI, Instead open \
    via browser manually.')
    args = parser.parse_args()
    module_path = plasmapy.__path__[0]
    calculator_path = os.path.dirname(module_path)
    computed_calculator_path = os.path.join(calculator_path,'docs','notebooks','plasma_calculator.ipynb')
    os.system(f'voila {computed_calculator_path}')
