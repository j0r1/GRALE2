.. _graleeditor:

GRALE Editor
============

TODO: for now these are just some very basic notes, mainly for myself. Needs
much cleaning up.

Command line
------------

Usage::

    grale_editor [arguments]

Is:
 - ``--nocheck`` : when exiting, don't ask to save changes if there are any

Starts with:
 - ``--zoom:`` : sets zoom level
 - ``--imgdataname:`` : sets latest images data file name, can be useful to set export name
 - ``--imgplane:`` : pre-loads the following file to be used in the back-projection related
   tools. This file should be a pickled ImagePlane or MultiImagePlane

Argument ends with 

 - .json: loads scene if first argument, adds to existing scene if not the 
   first argument
 - .png or .jpg: adds image layer with this image
 - .fits: adds fits layer for this file
 - .imgdat or .imgdata: adds the information for the images data file.
   By default, all images are added into a single image layer. By
   specifying a number followed by a comma and then the file name, that
   specific image is loaded into an image layer. If the number is negative,
   each image in the file is added as a separate image layer.
 
If the argument is just ``.``, an empty points layer is added.

Keyboard and mouse controls
---------------------------

Single click:
   - mouse on point: toggles point selection
   - adds point if points layer is active, match point if FITS/image 
     layer is active
   - clears selected points unless ctrl is pressed

If points layer active, and left mouse press and move:
   - if not on point: draws line along which points
     will be added when mouse is no longer held,
     triangles will be added as well. Interior points
     will be included in triangulation
   - if on point, that point will be moved. If control
     is pressed as well, other selected points will
     be moved as well

Right-click and move: select points/match points

Right-click:
 - on point: bring up dialog allowing you to set specific
   coordinates, time delay info and point group name
 - on FITS image: bring up dialog allowing you to change center and
   min/max value for the brightness scale

Double click point/match point: allows to edit point group name or match 
point name

Double click on one of a set of selected points: create triangulation

Double click elsewhere:
 - If FITS layer active: center on that point
 - If image layer active: match image to visible FITS layer
 - If points layer active: start contour finder based on what's visible
   around the clicked position. 

Control-delete: remove selected triangles

Shift-delete: remove selected points, and affected triangles

Control-C: copy

Control-X: cut

Control-Z: undo

Control-L: detect contour levels around central position (cfr double click)

Control-E: export visible images to images data file

Shift-Control-Z: redo

Double click on layer in list widget: make that layer active

Right click on layer in list widget: popup menu with options

Just 'c': center on selected points

Number 0-9:
 - No 'control' or 'alt': set zoom to 2^number
 - With 'control, but no 'alt': center point layer
