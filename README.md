# los_simulation
A repository to model line-of-site power from one or more transmitters and receivers, includes antenna directivity. 


# Installation 
```shell
pip install git+https://github.com/sandersa-nist/los_simulation.git
```
![image](./los_simulation/resources/example_scenario.gif)


# Example
```python
from los_simulation import *
create_scenario_1()
```
# Typical Workflow
 1. Create one or more RxNodes with a set location, direction and antenna pattern
 2. Create one or more TxNodes with a set location, direction, antenna pattern and power
 3. Use node_to_node_power to calculate the friis equation for each node of interest.
    + power_list_rx1 = np.array(list(map(lambda x: node_to_node_power(rx1,x,wavelength=C/frequency),txs)))
    + total_power_rx1 = 10*np.log10(np.sum(10**(power_list_rx1/10)))

    
## Jupyter Example
[Example](./examples/los_simulation_example.ipynb)

# API Documentation
[Start](./documentation/index.html)
