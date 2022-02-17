# dumbbellRheology

## Introduction
The behaviour of viscoelastic fluids are frequently simulated using a collection of elastic dumbbells suspended in a solvent. The stretching and relaxation of the dumbbells approximate the elastic behaviour of long polymeric molecules, making the polymeric solution exhibit both viscous and elastic properties. This repository contains scripts to simulate the viscoelastic behaviour of a suspension of dumbbells, which represent a polymeric liquid, in simple steady flows like shear and elongational flows. The nondimensionalized conformation tensor equations for such a solution, shown in [Bird, Robert Byron, et al. Dynamics of polymeric liquids, volume 2: Kinetic theory. Wiley, 1987](https://orbit.dtu.dk/en/publications/dynamics-of-polymeric-liquids-volume-2-kinetic-theory-2nd-edition) are solved by decomposing the symmetric conformation tensor into six coupled algebraic equations (six ODEs for unsteady state).

## Prequisites
Requires `scipy`

## Use

To test the behaviour of a dumbbell suspension with finite extensibility (b) = 50 and a relaxation time (t_rel) of 0.1, an instance of suspension should be created as :
```python
A=suspension(0.1,50)
```
The `eq_state()` method returns the conformation of the dumbbells at equilibrium

```python
A.eq_state()
```
The `shear(s)` method applies a specified steady shear of `s` on the suspension

```python
A.shear(100)
```
The `extension(e)` method applies a specified steady biaxial(for negative `e`) or uniaxial(for positive `e`) elongation on the suspension.

```python
A.shear(100)
```

The `QQ()` method returns the current conformation of the dumbbells in the suspension

```python
A.QQ()
```
