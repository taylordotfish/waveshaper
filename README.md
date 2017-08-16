Fish Waveshaper
===============

Fish Waveshaper is a powerful LV2 waveshaper plugin.


Usage
-----

Toggle the “Show plot” switch to display a graph of the current waveshaping
function. *You must have [gnuplot] installed for the graph to appear.* The
x-axis represents input amplitude and the y-axis represents output amplitude.
(Specifically, the axes represent the absolute value of the amplitude—the
waveshaping function is applied symmetrically to positive and negative
amplitudes.)

[gnuplot]: http://gnuplot.info/

You can then change the X and Y values for points 1 through 6 to add points to
the graph. The Y value of the rightmost point on the graph can be changed with
the “End point Y” control.


### Curves

The “curve”, “tension”, and “tangent” options for each point affect the
shape of the curve between that point and the next point to the right. The
“curve” option affects what type of curve is used:

* **Linear:** A straight line between the two points.
* **Tension:** A parabolic curve affected by the “tension” option. It is
  always increasing or decreasing; it does not change direction.
* **Tangents:** A curve with the slope of the leftmost part determined by
  “tan.&nbsp;1” and the rightmost part determined by “tan.&nbsp;2”.

The “Origin” tension, tangent, and curve settings affect the curve between the
point in the bottom left corner (representing an input and output amplitude of
0) and the next point to the right.


### Asymmetric mode

In asymmetric mode (enabled with the “Asymmetric” switch), the origin is no
longer in the bottom left corner, but is rather fixed in the center of the
graph. It still represents an input and output amplitude of 0, but now X values
to the left of the point represent negative input amplitude, and Y values below
the point represent negative output amplitude. (Values on the graph are scaled
to range from -1 to 1.) The top and bottom parts of the wave can now be shaped
differently.

With the origin now in the center of the graph, the “Asym.” options affect the
point in the lower-left corner. “Asym. Y” affects the Y value of this point,
and the “Asym.” tension, tangent, and curve settings affect the curve between
this point and the next point to the right.


### Other options

“Monitor levels” enables the use of the peak meters, and allows you to see
input level, output level, and whether or not the input is being clipped. (The
input is clipped to 0dBFS, which corresponds to an amplitude of 1 or -1.)

“Oversample” enables oversampling by the specified factor, which reduces
aliasing. This can sometimes amplify the output signal (you may want to use the
“Output gain” control in this case).

“Input gain” controls the amplification of the input signal, and “Output gain”
controls the amplification of the output signal.


Dependencies
------------

* LV2 development files
* gnuplot
* GCC
* GNU Make

On Debian GNU/Linux (and many derivatives), these can be installed by running
the following command as root:

```
apt-get install lv2-dev gnuplot gcc make
```


Installation
------------

Run the following commands (you will need to have [Git] installed):

```
git clone https://github.com/taylordotfish/waveshaper ~/.lv2/fish-waveshaper.lv2/
cd ~/.lv2/fish-waveshaper.lv2/
make
```

[Git]: https://git-scm.com/


License
-------

Fish Waveshaper is licensed under the GNU General Public License, version 3 or
any later version. See [LICENSE].

This README file has been released to the public domain using [CC0].

[LICENSE]: LICENSE
[CC0]: https://creativecommons.org/publicdomain/zero/1.0/
