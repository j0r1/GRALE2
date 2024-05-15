.. _notebooks:

Example Jupyter notebooks
=========================

 * Example plots of time-delay surfaces: `tdsurfaces.ipynb <_static/tdsurfaces.ipynb>`_

 * Using a previously calculated deflection angle grid to create a 
   :class:`DeflectionGridLens <grale.lenses.DeflectionGridLens>`: 
   `deflectiongridlens.ipynb <_static/deflectiongridlens.ipynb>`_

 * The `LensPerfect <http://adsabs.harvard.edu/abs/2008ApJ...681..814C>`_
   method, using :class:`Wendland lenses <grale.lenses.MultipleWendlandLens>`
   to obtain a mass distribution if deflection angles are known
   at various locations: `wendland.ipynb <_static/wendland.ipynb>`_

 * Illustrations of using a grid function, as well as of the subdivision grid: 
   `gridtests.ipynb <_static/gridtests.ipynb>`_

 * Illustration of how to approximate a mass distribution by multiple Plummer
   basis functions: `fittest.ipynb <_static/fittest.ipynb>`_

 * Creating a circularly symmetric gravitational lens based on the
   specification of it's profile: `profilelens.ipynb <_static/profilelens.ipynb>`_

 * Illustration of the mass sheet degeneracy in its simplest guise: 
   `msdexample.ipynb <_static/msdexample.ipynb>`_

 * An example that shows explicitly how you can modify an existing
   lens model to influence only the time delay: `timedelayadjust.ipynb <_static/timedelayadjust.ipynb>`_

 * Example notebook that illustrates how you can create a degenerate lens
   that moves the source, by combining two mass disk degeneracies:
   `massdisk_movesource_smooth.ipynb <_static/massdisk_movesource_smooth.ipynb>`_

 * An example with multiple lens planes; also shows the effect as
   the source redshift increases: `multilensplane.ipynb <_static/multilensplane.ipynb>`_

 * Multi-lensplane example from 
   `Compound lensing: Einstein Zig-Zags and high multiplicity lensed images <http://adsabs.harvard.edu/abs/2016MNRAS.456.2210C>`_

   * Situation in Figure 7: `multisistests.ipynb <_static/multisistests.ipynb>`_
   * Situation in Figure 2: `multisistests2.ipynb <_static/multisistests2.ipynb>`_
   * Situation in Figure 6 (second row): `multisistests3.ipynb <_static/multisistests3.ipynb>`_

 
 * Illustration of loading `LensTool <https://projets.lam.fr/projects/lenstool/wiki>`_
   models: `lenstooltest.ipynb <_static/lenstooltest.ipynb>`_

 * An example of generating fake weak lensing measurements: 
   `generatefakewldata.ipynb <_static/generatefakewldata.ipynb>`_

 * Recreates some plots of the article `A generalization of the mass-sheet degeneracy 
   producing ring-like artefacts in the lens mass distribution <https://ui.adsabs.harvard.edu/abs/2008MNRAS.386..307L/abstract>`_,
   describing a generalization of the mass sheet degeneracy that works with sources
   at different redshifts: `scaledegen.ipynb <_static/scaledegen.ipynb>`_

 * The previous example scales two sources at different redshifts with the
   same scale factor, in the article `Lensing degeneracies and mass substructure <https://ui.adsabs.harvard.edu/abs/2012MNRAS.425.1772L/abstract>`_
   it was shown how different scale factors can be used instead. Some results
   can be found here: `scaledegen2012.ipynb <_static/scaledegen2012.ipynb>`_

 * These notebooks show how the IrtyshI and IrtyshII lens models, used in
   `Free-form grale lens inversion of galaxy clusters with up to 1000 multiple images <https://ui.adsabs.harvard.edu/abs/2020MNRAS.494.3998G/abstract>`_
   and created using the `gravlens/lensmodel software <https://www.physics.rutgers.edu/~keeton/gravlens/2012WS/>`_
   (see also `Keeton 2001 <https://ui.adsabs.harvard.edu/abs/2001astro.ph..2341K/abstract>`_),
   can be convert to lens models for Grale: `irtyshI.ipynb <_static/irtyshI.ipynb>`_ and
   `irtyshII.ipynb <_static/irtyshII.ipynb>`_

 * Importing deflection fields from `various models of the Abell 370 <https://archive.stsci.edu/pub/hlsp/frontier/abell370/models/>`_:

   * `a370cats.ipynb <_static/a370cats.ipynb>`_
   * `a370diego.ipynb <_static/a370diego.ipynb>`_
   * `a370glafic.ipynb <_static/a370glafic.ipynb>`_
   * `a370keeton.ipynb <_static/a370keeton.ipynb>`_
   * `a370sharon.ipynb <_static/a370sharon.ipynb>`_
   * `a370williams.ipynb <_static/a370williams.ipynb>`_

 * A modification of `msdexample.ipynb <_static/msdexample.ipynb>`_ above, to illustrate the
   generation of equivalent lens models by extrapolating the lens potential: `msdexample-equivlenstests.ipynb <_static/msdexample-equivlenstests.ipynb>`_

 * Using the code from the lens potential extrapolation to obtain lenses with different
   MSD-like effects for different sources: `potentialextrap_multisheet.ipynb <_static/potentialextrap_multisheet.ipynb>`_

 * This example uses the A2744 data from `Weak gravitational lensing measurements of Abell 2744 using JWST and shear measurement algorithm pyRRG-JWST <https://ui.adsabs.harvard.edu/abs/2024MNRAS.529..802H/abstract>`_
   to recreate (more or less) their Fig. 6 plot that estimates the mass density from the weak
   lensing measurements: `a2744-wldatatest.ipynb <_static/a2744-wldatatest.ipynb>`_

