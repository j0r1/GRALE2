Lens inversion examples
-----------------------

These examples are based on the ones for the [original](http://research.edm.uhasselt.be/jori/page/Physics/GraleV1.html)
version of Grale/Graleshell. The way the Python scripts are built
is fairly straightforward. Of course, because of the use of Python now,
things can be automated in many more interesting ways than is done
in these scripts, but for the purpose of illustrating the main commands
it is probably not a bad idea to keep things simple.

In each example, the procedure is very similar:

 - set up the inversion (cosmological model, input files, etc)
 - start from a uniform grid for a first inversion
 - perform the inversion
 - based on the mass distribution that's retrieved, create a
   subdivision grid
 - this inversion/subdivision combination is repeated a number
   of times, with an increasing amount of grid cells to allow
   for a more detailed reconstruction

The lens models that are reconstructed along the way are saved to
files, and they can be visualized in the usual way, as explained
in the [tutorial](http://research.edm.uhasselt.be/~jori/grale2/tutorial.html)
for example.

From such a set of lens models that's generated using this subdivision
procedure, you should then choose the one that you think is the best
solution. This is **not** done automatically in these scripts, one
reason is that this would complicate the scripts themselves somewhat
and I like to keep them simple to illustrate the core procedure.
Another reason is that, while determining the best solution of the set
is trivial in case there's just one fitness criterion, it gets more
complicated when using more that one fitness measure. Sometimes eye-balling
the fitness values of the reconstructions can help to discard some
runs that got stuck in some local minimum that clearly does not represent
something even close to a solution of the inversion.

After choosing the best model from a set of these models, typically
you'd repeat the whole procedure (so execute the `inversion.py` script
a number of times and select a model), to get an idea about the average
solution and the dispersion for example. Again, these results can be
visualized in the usual way.
