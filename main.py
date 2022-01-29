import json
from entities.simulation import Simulation

with open('./configs/system_configs.json', 'r') as f:
    sys_configs = json.load(f)

with open('./configs/process_configs.json', 'r') as f:
    process_configs = json.load(f)

with open('input.json', 'r') as f:
    input = json.load(f)


simulation = Simulation(input,sys_configs,process_configs)
simulation.run_simulation()

print('fim')