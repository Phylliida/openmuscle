# openmuscle
An open source implementation of the Hill-type muscle model without any dependencies.

Example usage (in python)

```
import openmuscle
activation = 0.3 # Activation may be between 0.001 and 1 (anything else is clamped to these values)
muscle = openmuscle.Muscle(activation)

# Simulate the muscle for 0.1 seconds 
# This will output the resulting force and update the muscle
print muscle.step(activation, 0.1)
```
