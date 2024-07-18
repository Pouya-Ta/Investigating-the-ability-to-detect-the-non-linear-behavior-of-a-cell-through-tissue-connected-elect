# Investigating the Ability to Detect Non-Linear Behavior of a Cell through Tissue-Connected Electrodes

`main.py` is a Python script designed to investigate the non-linear behavior of a cell by simulating the Hodgkin-Huxley (HH) model of a neuron. The script explores the effects of various parameters on membrane potential and the gating variables of potassium (K) and sodium (Na) channels. It also calculates equilibrium points of the HH model for different current injection values.

## Features

- **Hodgkin-Huxley Model Simulation**: Implements the HH model to simulate the neuron's behavior under different conditions.
- **Parameter Variation**: Allows for the variation of injected current and leakage voltage to observe their effects on membrane potential and gating probabilities.
- **Equilibrium Points Calculation**: Computes the equilibrium points of the HH model for different current injection values.
- **Visualization**: Provides plots of membrane potential and gating probabilities over time for various conditions.

## How It Works

1. **Initialization**: The script initializes constants related to the HH model, such as conductances, reversal potentials, membrane capacitance, and temperature coefficient.
2. **Rate Functions**: Defines the opening (alpha) and closing (beta) rate functions for the K and Na channels.
3. **Simulation**: 
   - Simulates the HH model for a range of injected current values, plotting the membrane potential and gating probabilities over time.
   - Examines the effect of varying leakage voltage on membrane potential.
   - Calculates and plots equilibrium points for different current injection values.
4. **Visualization**: Uses `matplotlib` to generate plots that illustrate the neuron's response under various conditions.

## Usage

1. Ensure you have the required libraries installed:
    ```bash
    pip install numpy matplotlib scipy sympy
    ```
2. Run the script:
    ```bash
    python main.py
    ```
3. Follow the output plots to analyze the neuron's behavior under different conditions.

## Example Output

### Membrane Potential and Gating Probabilities
For each current injection value, the script generates two plots:
- **Membrane potential vs. time**
- **Gating probabilities (n, m, h) vs. time**

![Example Output](example_output.png)

### Effect of Varying Leakage Voltage
The script also generates a plot showing the effect of different leakage voltage values on membrane potential.

![Leakage Voltage Effect](leakage_voltage_effect.png)

### Equilibrium Points
Equilibrium points for the HH model are computed and printed for a range of current injection values, with corresponding plots showing membrane potential over time.

![Equilibrium Points](equilibrium_points.png)

## Code Overview

### Constants and Parameters
```python
# Constants and parameters initialization
gK = 36  # unit: mS/cm^2
vK = -72  # unit: mV
gNa = 120  # unit: mS/cm^2
vNa = 55  # unit: mV
gL = 0.3  # unit: mS/cm^2
vL = -49.4  # unit: mV
C = 1  # unit: muF/cm^2
phi = 1
Q = 3  # Temperature coefficient-unitless
```

### Rate Functions
```python
# Rate functions for K and Na channels
an = lambda V: (-0.01*(V+50))/(np.exp(-(V+50)/10)-1)
bn = lambda V: 0.125*np.exp(-(V+60)/80)
am = lambda V: (-0.1*(V+35))/(np.exp(-(V+35)/10)-1)
bm = lambda V: 4*np.exp(-(V+60)/18)
ah = lambda V: 0.07*np.exp(-(V+60)/20)
bh = lambda V: 1/(np.exp(-(V+30)/10)+1)
```

### Simulation and Visualization
```python
# Simulation and plotting for various current injection values
for I in range(0, 11, 1):
    def func(p, t):
        return [(1/C)*(I-(gK*p[1]**4*(p[0] - vK)) - (gNa*p[2]**3*p[3]*(p[0]-vNa))-(gL*(p[0]-vL))),
                (an(p[0])*(1-p[1])) - (bn(p[0])*p[1]),
                (am(p[0])*(1-p[2])) - (bm(p[0])*p[2]),
                (ah(p[0])*(1-p[3])) - (bh(p[0])*p[3])]
    Parameters = odeint(func, y0=[-60, 0.317, 0.0529, 0.596], t=tspan)
    plt.figure()
    plt.plot(tspan, Parameters[:,0], color='red')
    plt.xlabel('time')
    plt.ylabel('membrane Voltage (mV)')
    plt.title(f'Membrane potential/time, I = {I} μA/cm^2')
    plt.show()
    plt.figure()
    plt.plot(tspan, Parameters[:,1], label='n')
    plt.plot(tspan, Parameters[:,2], label='m')
    plt.plot(tspan, Parameters[:,3], label='h')
    plt.legend()
    plt.xlabel('t')
    plt.title(f'n, m & h/t, I = {I} μA/cm^2')
    plt.show()
```
