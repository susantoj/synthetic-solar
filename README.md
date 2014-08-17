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

References
----------

+ R. Aguiar and M. Collares-Pereira, ["A simple procedure for the generation of sequences of daily radiation values using Markov transition matrices"](http://www.sciencedirect.com/science/article/pii/0038092X88900497), Solar Energy, vol. 40, 269-279, 1988
+ R. Aguiar and M. Collares-Pereira, ["TAG: A time-dependent, autoregressive Gaussian model for generating synthetic hourly radiation"](http://www.sciencedirect.com/science/article/pii/0038092X9290068L), Solar Energy, vol. 49, 167-174, 1992
+ A. Luque, S. Hegedus, [“Handbook of Photovoltaic Science and Engineering”](http://www.wiley.com/WileyCDA/WileyTitle/productCd-0470721693.html), Wiley, 2003

Credits
-------

+ **Julius Susanto** - http://github.com/susantoj

License
-------

This program is free software: you can redistribute it and/or modify
it under the terms of the [GNU General Public License](http://www.gnu.org/copyleft/gpl.html) as published
by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.