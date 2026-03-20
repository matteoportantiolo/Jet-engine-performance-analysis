# Jet Engine Performance Analysis (Jumo 004 vs BMW 003)

## Overview
This project presents a detailed analysis and comparison of two early turbojet engines: the Jumo 004 and the BMW 003.

The study focuses on thermodynamic cycle analysis, compressor design, and performance enhancement strategies, providing a comprehensive evaluation of propulsion system behavior across different flight conditions.

## Objectives
- Analyze and compare turbojet engine performance  
- Compute thermodynamic cycle properties  
- Design and analyze axial compressors  
- Evaluate thrust augmentation strategies  

## Methods

### Thermodynamic Cycle Analysis
- Full cycle computation for both engines  
- Evaluation across a range of flight velocities  
- Validation against literature data  

At a fixed altitude of 30,000 ft, atmospheric conditions were modeled and used to compute engine performance :contentReference[oaicite:0]{index=0}  

### Compressor Design
- Multi-stage axial compressor analysis  
- Velocity triangles and stage-by-stage modeling  
- Free vortex theory for blade twist  

The compressor geometry and performance were derived using engineering assumptions such as constant axial velocity and repeated stage modeling :contentReference[oaicite:1]{index=1}  

### Performance Comparison
- Thrust  
- Specific fuel consumption (TSFC)  
- Thermal and propulsive efficiency  

Results show similar performance between the two engines, with slightly higher thrust for the Jumo 004 due to higher compression ratio :contentReference[oaicite:2]{index=2}  

### Thrust Augmentation Strategies
Two different approaches were analyzed:

- **Afterburner (Jumo 004)**  
- **Liquid Rocket Engine (BMW 003)**  

The LRE system was designed including:
- Combustion chamber sizing  
- Injection system  
- Feeding system (turbopump vs pressurized tank)  

## Results

Key findings:

- Jumo 004 provides slightly higher thrust (~2.5% at design condition) :contentReference[oaicite:3]{index=3}  
- BMW 003 shows comparable efficiency with lower compression ratio  
- Compressor design strongly affects overall performance  
- Thrust augmentation significantly increases performance at the cost of efficiency  

## Implementation

The project is implemented in MATLAB and includes:

- Thermodynamic cycle solvers  
- Compressor stage analysis  
- Performance comparison tools  
- Propulsion system sizing  

## Key Concepts
- Turbojet engines  
- Thermodynamic cycles  
- Axial compressors  
- Propulsion systems  
- Afterburner and rocket augmentation  

## Author
Matteo Portantiolo  
MSc Space Engineering – GNC
