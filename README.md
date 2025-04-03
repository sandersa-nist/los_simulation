# los_simulation
A repository to model line-of-site power from one or more transmitters and receivers, includes antenna directivity. 
Based on the [friis](https://en.wikipedia.org/wiki/Friis_transmission_equation) equation. 

# Installation 
```shell
pip install git+https://github.com/sandersa-nist/los_simulation.git
```
# Example Scenario
```python
from los_simulation import *
create_scenario_1()
```
# Typical Workflow
 1. Create one or more RxNodes with a set location, direction and antenna pattern
 2. Create one or more TxNodes with a set location, direction, antenna pattern and power
 3. Use node_to_node_power to calculate the friis equation for each node of interest.
 
 ```python
from los_simulation import *

frequency =3e9
rx = RxNode(location=[0,0],direction=[0,1],antenna_pattern=simple_directional_gain)
txs = []
number_transmitters =11
tx_distance = 1000
tx_angles = np.linspace(-np.pi,np.pi-2*np.pi/number_transmitters,number_transmitters)
for i,angle in enumerate(tx_angles):
    location = [tx_distance*np.cos(angle),tx_distance*np.sin(angle)]
    new_tx = TxNode(location=location,direction = [0,1],antenna_pattern = omni, power = -10, id=f"tx_{i}")
    txs.append(new_tx)
power_list_rx = np.array(list(map(lambda x: node_to_node_power(rx,x,wavelength=C/frequency),txs)))
total_power_rx = 10*np.log10(np.sum(10**(power_list_rx/10)))
```
# [Jupyter Example](./examples/los_simulation_example.ipynb)
# [API Documentation](./documentation/index.html)

![image](./los_simulation/resources/example_scenario.gif)