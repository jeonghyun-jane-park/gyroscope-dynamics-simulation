# Gyroscope Dynamics Simulation & Experimental Validation

A nonlinear rigid-body dynamics model of a laboratory gyroscope, implemented in MATLAB and compared against measured angular velocity and acceleration from a body-mounted sensor.

This project was originally completed as a team assignment for **MCEN90038 Dynamics at the University of Melbourne (May 2026)** by **Jeonghyun Park, Sam Seaberry, and Katherine Tao**. This repository reorganizes the original work into a cleaner technical portfolio format while preserving team attribution.

## Project overview

The model captures three coupled rotational degrees of freedom:

- **Precession** about the vertical axis, `beta`
- **Tilt / nutation** of the main bar, `theta`
- **Rotor spin** about the rotor axis, `alpha`

The MATLAB workflow:

1. defines successive body-frame rotation matrices;
2. constructs inertia tensors for cylindrical components with central rod cutouts;
3. transfers inertia tensors to the pivot with the parallel-axis theorem;
4. applies Newton-Euler angular-momentum balance in rotating frames;
5. symbolically solves for the three angular accelerations;
6. integrates the nonlinear equations of motion with `ode45`;
7. reconstructs angular velocity and linear acceleration at the sensor location; and
8. compares the simulated response with experimental measurements.

## Key engineering concepts

- 3D rigid-body kinematics
- Rotation matrices and body-fixed frames
- Angular momentum and the transport theorem
- Newton-Euler dynamics
- Parallel-axis theorem
- Gyroscopic precession and nutation
- Symbolic equation generation in MATLAB
- Numerical integration with `ode45`
- Sensor-frame acceleration reconstruction
- Experimental model validation

## Repository structure

```text
gyroscope-dynamics-simulation/
├── README.md
├── src/
│   └── gyroscope_simulation.m
├── report/
│   ├── gyroscope_dynamics_report.tex
│   ├── gyroscope_dynamics_report.pdf
│   └── figures/
├── data/
│   └── README.md
└── media/
    └── README.md
```

## Model formulation

The generalized coordinates are

```text
q = [beta, theta, alpha]^T
```

and the symbolic equations are assembled into the nonlinear form

```text
A(q) q_ddot = b(q, q_dot)
```

The code solves this system symbolically for the angular accelerations, converts the resulting expressions to numerical MATLAB functions, and integrates the six-state system over time.

For the bar/sensor frame, the angular velocity is reconstructed as

```text
omega_2 = [theta_dot,
           beta_dot sin(theta),
           beta_dot cos(theta)]^T
```

and the sensor acceleration includes tangential, centripetal, and gravitational contributions.

## Initial conditions used in the original experiment

- Rotor speed: approximately **500 rpm** (`52.36 rad/s`)
- Initial bar tilt: **15 deg** under the model coordinate convention
- Other initial angular rates: **0 rad/s**
- Simulation duration: **90 s**

The model assumes a rigid, frictionless system with a massless connecting rod. These assumptions explain part of the long-time mismatch between simulation and experiment.

## Results

The simulation reproduces the main qualitative features of the physical gyroscope:

- the same precession direction relative to rotor spin;
- nutation-like oscillation following release;
- similar short-time angular-velocity amplitudes and periodic structure; and
- broadly similar sensor-acceleration behavior over short windows.

The largest long-duration discrepancies are consistent with unmodeled friction/drag, uncertainty in the initial rotor speed, simplified geometry, and possible sensor/alignment effects. The revised report discusses these limitations in more detail.

## Running the MATLAB model

### Requirements

- MATLAB
- Symbolic Math Toolbox

Run:

```matlab
run('src/gyroscope_simulation.m')
```

The script can generate the numerical simulation and animation without the experimental CSV. To reproduce the measured-vs-calculated plots, place the original sensor export at:

```text
data/ExportedData - Run 2.csv
```

The expected column names are documented in [`data/README.md`](data/README.md).

> **Data note:** The raw experimental CSV was not included in the files used to build this public-facing repository package, so no synthetic or reconstructed dataset has been added.

## Technical report

The polished LaTeX report is available in:

- [`report/gyroscope_dynamics_report.pdf`](report/gyroscope_dynamics_report.pdf)
- [`report/gyroscope_dynamics_report.tex`](report/gyroscope_dynamics_report.tex)

The revised version removes student ID numbers, improves the mathematical presentation, adds a parameter table and reproducibility note, and organizes the work around modeling, implementation, validation, and limitations.

## Possible next steps

A stronger validation study could add:

- friction or drag terms identified from measured rotor-speed decay;
- an exact mass-weighted center of mass for the rotor/pulley assembly;
- parameter estimation for uncertain initial conditions;
- RMSE and phase-error metrics for each sensor channel; and
- automated generation of publication-quality plots directly from the raw CSV.

## Attribution

This was a **team project**. Public portfolio use should preserve the names of all original contributors unless the team agrees on another attribution format. Individual contribution statements can be added later if each contributor's role is known and accurately documented.
