Synthetic Solar
===============

Synthetic Solar is a Python module for generating synthetic sequences of hourly
solar irradiances.

Example
-------

```python 
import synth_solar
import matplotlib.pyplot as plt

Ktm = [0.49, 0.48, 0.50, 0.51, 0.51, 0.54, 0.54, 0.57, 0.56, 0.54, 0.51, 0.51]
lat = -34.8
G0 = synth_solar.Aguiar_hourly_G0(Ktm, lat)

plt.plot(G0)
plt.show()
```

![screenshot of sample output](/example.png?raw=true)

Credits
-------

+ **Julius Susanto** - http://github.com/susantoj

License
-------

This program is free software: you can redistribute it and/or modify
it under the terms of the [GNU General Public License](http://www.gnu.org/copyleft/gpl.html) as published
by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.