# Simulation of Neuronal Behavior Using Hodgkin-Huxley Model

This repository contains both Python (`main.py`) and MATLAB scripts for simulating the non-linear behavior of a neuron based on the Hodgkin-Huxley (HH) model. The scripts explore the dynamics of membrane potential and the conductances of potassium (K) and sodium (Na) channels under various conditions, such as applied current and voltage clamping.

## Features

- **Hodgkin-Huxley Model Simulation**: Implements the HH model to simulate the dynamics of a neuron, focusing on the conductances of Na and K channels and their impact on membrane potential.
- **Voltage-Clamp Experiments**: Simulates voltage-clamp experiments to observe how potassium and sodium conductances evolve over time for a given membrane potential.
- **Action Potential Generation**: Models the generation of action potentials by applying an external current to the neuron and visualizing the resulting membrane potential and channel conductances.
- **Visualization**: Provides detailed plots of membrane potential, sodium and potassium conductances, and gating variables over time.

## MATLAB Code Overview

### Usage

1. **Run the MATLAB script**:
   - Open the MATLAB script (`voltage_clamp_HH.m`) in MATLAB.
   - Execute the script by pressing the "Run" button or typing `run('voltage_clamp_HH.m')` in the command window.
   
2. **Output**:
   - The output plots will provide insights into the dynamics of membrane potential, potassium and sodium conductances, and gating variables.

### Example Outputs

- **Voltage-Clamp Simulation**:
  - **Potassium Conductance**:
    ![MATLAB Potassium Conductance](matlab_potassium_conductance.png)
  - **Sodium Conductance**:
    ![MATLAB Sodium Conductance](matlab_sodium_conductance.png)

- **Action Potential Simulation**:
  - **Membrane Potential**:
    ![MATLAB Action Potential](matlab_action_potential.png)
  - **Sodium and Potassium Conductances**:
    ![MATLAB Conductances](matlab_conductances.png)

### Code Overview

#### Constants and Parameters

```matlab
% Constants and parameters initialization
gK = 36;    % Potassium conductance (mS/cm^2)
gNa = 120;  % Sodium conductance (mS/cm^2)
gL = 0.3;   % Leak conductance (mS/cm^2)
vK = -12;   % Potassium reversal potential (mV)
vNa = 115;  % Sodium reversal potential (mV)
vL = 10.6;  % Leak reversal potential (mV)
C_m = 1e-6; % Membrane capacitance (F)
```

#### Differential Equations

```matlab
function dXdt = dXdT_HH(t, X, I_app)
    v = X(1); m = X(2); n = X(3); h = X(4);
    % Differential equations governing the system
    ...
    dXdt = [dv_dt; dm_dt; dn_dt; dh_dt];
end
```

#### Simulation and Visualization

```matlab
[t, X] = ode15s(@(t, X) dXdT_HH(t, X, I_app), [0, 30], [0, 0, 0, 0]);
figure;
plot(t, X(:, 1), 'k-', 'LineWidth', 1.5);
xlabel('t (ms)');
ylabel('v (mV)');
title('Membrane Potential Over Time');
```

## Python Code Overview

### Usage

1. **Ensure the required Python libraries are installed**:

   ```bash
   pip install numpy matplotlib scipy
   ```

2. **Run the script**:

   ```bash
   python main.py
   ```

3. **Output**:
   - The output plots will provide insights into the dynamics of membrane potential, potassium and sodium conductances, and gating variables.

### Example Outputs

- **Voltage-Clamp Simulation**:
  - **Potassium Conductance**:
    ![Python Potassium Conductance](python_potassium_conductance.png)
  - **Sodium Conductance**:
    ![Python Sodium Conductance](python_sodium_conductance.png)

- **Action Potential Simulation**:
  - **Membrane Potential**:
    ![Python Action Potential](python_action_potential.png)
  - **Sodium and Potassium Conductances**:
    ![Python Conductances](python_conductances.png)

### Code Overview

#### Constants and Parameters

```python
# Constants and parameters initialization
gK = 36    # Potassium conductance (mS/cm^2)
gNa = 120  # Sodium conductance (mS/cm^2)
gL = 0.3   # Leak conductance (mS/cm^2)
vK = -12   # Potassium reversal potential (mV)
vNa = 115  # Sodium reversal potential (mV)
vL = 10.6  # Leak reversal potential (mV)
C_m = 1e-6 # Membrane capacitance (F)
```

#### Differential Equations

```python
def dXdT_HH(t, x, I_app):
    v, m, n, h = x
    # Differential equations governing the system
    ...
    return [dv_dt, dm_dt, dn_dt, dh_dt]
```

#### Simulation and Visualization

```python
sol = solve_ivp(lambda t, x: dXdT_HH(t, x, I_app)[0], [0, 30], [0, 0, 0, 0], method='BDF')
plt.plot(t, x[:, 0], 'k-', linewidth=1.5)
plt.xlabel('t (ms)')
plt.ylabel('v (mV)')
plt.title('Membrane Potential Over Time')
plt.show()
```

## Comparison: MATLAB vs. Python

### Pros of MATLAB

1. **Ease of Use**: MATLAB is designed for mathematical and engineering computations, offering built-in functions for numerical integration (`ode23s`, `ode15s`) that are simple to use with minimal setup.
2. **Performance**: MATLAB’s numerical solvers are highly optimized for handling stiff ODEs, making it particularly efficient for solving complex biological models like Hodgkin-Huxley.
3. **Visualization**: MATLAB has powerful plotting tools that are highly customizable and widely used in scientific research, which is beneficial for creating publication-quality figures.

### Cons of MATLAB

1. **Cost**: MATLAB is a proprietary software, requiring a paid license, which can be a significant barrier for many users, especially in academia or small research groups.
2. **Less Flexibility**: While MATLAB is powerful, it is less versatile than Python for general-purpose programming and integrating with other technologies, such as web applications or machine learning frameworks.

### Pros of Python

1. **Open Source**: Python is free and open-source, making it accessible to a wide audience. It is supported by a large community, ensuring continuous improvements and extensive documentation.
2. **Flexibility**: Python is a general-purpose programming language that can be easily integrated with other technologies and tools, such as machine learning libraries (e.g., TensorFlow, PyTorch) and web frameworks (e.g., Django).
3. **Extensive Libraries**: Python’s ecosystem includes powerful libraries like `SciPy`, `NumPy`, and `Matplotlib`, which provide robust tools for scientific computing and data visualization.

### Cons of Python

1. **Performance**: Python can be slower than MATLAB for certain numerical computations, especially when dealing with large datasets or complex simulations. However, this can often be mitigated by using optimized libraries or parallel processing.
2. **Steeper Learning Curve for Scientific Computing**: While Python is a versatile language, setting up scientific computations might require more effort compared to MATLAB, especially for users unfamiliar with Python’s ecosystem.

### Conclusion

- **MATLAB** is ideal for users focused on mathematical modeling and simulations, particularly in engineering and biological sciences, due to its performance and ease of use in these domains.
- **Python** offers greater flexibility and is more accessible, making it a better choice for those looking to integrate simulations with other programming tasks or who require a free, open-source solution.

For this project, if you need to focus primarily on the simulation and visualization aspects, MATLAB might be more efficient. However, if you plan to extend the project to include machine learning, data analysis, or integration with other tools, Python would be the preferred option.
