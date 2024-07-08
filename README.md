# MoonHarbour

## Overview
MoonHarbour is an innovative project dedicated to executing rendezvous and docking maneuvers with the Lunar Gateway. This initiative plays a crucial role in advancing lunar exploration and establishing a sustainable human presence on the Moon.

## Lunar Gateway
The Lunar Gateway is an essential component of lunar exploration. It serves as a multi-purpose outpost orbiting the Moon, providing vital support for long-term human return to the lunar surface and as a staging point for deep space exploration. For more information, visit the [Lunar Gateway Wikipedia page](https://en.wikipedia.org/wiki/Lunar_Gateway).

## Rendezvous and Docking
Rendezvous and docking are two critical phases in spaceflight involving the approach and mechanical joining of two spacecraft. This complex process includes:
- **Rendezvous**: Maneuvering the spacecraft to approach each other in orbit.
- **Docking**: Physically connecting the spacecraft to transfer crew or cargo.

## Rendezvous Procedure Animation
![Rendezvous Procedure](Output/rendezvous.gif)

## Objectives of MoonHarbour
1. **Precise Maneuvering**: Develop and implement precise orbital maneuvering techniques for successful rendezvous.
2. **Innovative Docking Solutions**: Create and integrate advanced docking systems compatible with the Lunar Gateway.
3. **Safety and Reliability**: Ensure the highest standards of safety and reliability in all rendezvous and docking operations.
4. **Supporting Lunar Exploration**: Contribute to the broader goal of sustainable lunar exploration and long-term space habitation.

## Technical Details
This project builds on high-fidelity orbital dynamical modeling using the extended Battin-Giorgi equations to understand relative orbit motion accurately. Key contributions include:

### Dynamics and Control Strategies
- **Trajectory Control**: Implemented through feedback linearization algorithms to ensure collision-free maneuvers.
- **Attitude Control**: Ensured via feedback linearization and Lyapunov's method for stability and precise alignment.

### Actuation and Steering
- **Actuation Devices and Steering Laws**: Modeled to effectively implement control strategies.

### Trajectory Planning
- **Two-point and Three-point Methods**: Employed for safe trajectory planning and collision avoidance, including emergency maneuvers.

### Mission Phases
1. **Rendezvous (1500m to 100m)**: Focused on orbital control to correctly position the spacecraft.
2. **Final Approach (100m to 10m)**: Involves orbital control, attitude adjustments, and actuation for precise alignment.
3. **Drift and Mating (10m to 5m)**: Maintains precise attitude control and actuation for smooth docking or berthing.

### Attitude Control System
- Designed using a feedback law with quasi-global stability properties, ensuring correct spacecraft orientation for successful docking or berthing.

### Testing and Validation
- Extensive testing of docking and berthing scenarios under nominal and non-nominal conditions.
- Monte Carlo simulations demonstrate the effectiveness of control strategies for safe and precise maneuvers.

## Collaboration
MoonHarbour is open to collaboration with international space agencies, research institutions, and private sector partners.

## Contact
For more information or to discuss potential collaborations, please contact [me](mailto:leonardo.rxsso@gmail.com).
